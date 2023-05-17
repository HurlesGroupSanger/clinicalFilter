"""
Copyright (c) 2021 Genome Research Limited
Author: Ruth Eberhardt <re3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import logging

from utils.utils import add_compound_het_to_candidates
from utils.utils import add_single_var_to_candidates
from utils.utils import convert_genotype_to_gt


class AutosomalFilter(object):
    """
    Inheritance filters for SNVs/indels on autosomes
    """

    def __init__(self, Inheritancefilter, hgncid):
        self.family = Inheritancefilter.family
        self.parents = Inheritancefilter.parents
        self.variants_per_gene = Inheritancefilter.variants_per_gene
        self.candidate_variants = Inheritancefilter.candidate_variants
        self.inheritance_report = Inheritancefilter.inhreport
        self.gene = Inheritancefilter.genes[hgncid]
        self.hgncid = hgncid

    def autosomal_filter(self):
        if self.parents == "both":
            self.autosomal_both_parents()
        elif self.parents == "none":
            self.autosomal_no_parents()
        else:
            self.autosomal_single_parent()

    def autosomal_both_parents(self):
        """
        Screens variants in autosomes where there are both parents and adds
        candidates to candidate_variants
        """
        variants = self.variants_per_gene[self.hgncid]

        for v in variants.keys():
            if not variants[v]["child"].is_snv():
                continue
            mum_genotype = variants[v]["child"].get_mum_genotype()
            dad_genotype = variants[v]["child"].get_dad_genotype()
            mum_gt = convert_genotype_to_gt(mum_genotype)
            dad_gt = convert_genotype_to_gt(dad_genotype)
            mum_aff = self.family.mum.affected
            dad_aff = self.family.dad.affected

            if variants[v]["child"].genotype == "1":
                # heterozygous
                for inh in self.gene["mode"]:
                    if inh == "Biallelic":
                        self.biallelic_heterozygous_parents_filter(
                            v, variants[v]["child"], mum_gt, dad_gt, mum_aff, dad_aff
                        )
                    elif inh == "Monoallelic":
                        self.monoallelic_heterozygous_parents_filter(
                            v, variants[v]["child"], mum_gt, dad_gt, mum_aff, dad_aff
                        )
                    elif inh == "Mosaic":
                        self.mosaic_heterozygous_parents_filter(
                            v, variants[v]["child"], mum_gt, dad_gt, mum_aff, dad_aff
                        )
                    elif inh == "Imprinted":
                        self.imprinted_heterozygous_parents_filter(
                            v, variants[v]["child"], mum_gt, dad_gt, mum_aff, dad_aff
                        )
                    else:
                        logging.info(v + " unknown gene mode " + inh)
            elif variants[v]["child"].genotype == "2":
                # homozygous
                for inh in self.gene["mode"]:
                    if inh == "Biallelic":
                        self.biallelic_homozygous_parents_filter(
                            v, variants[v]["child"], mum_gt, dad_gt, mum_aff, dad_aff
                        )
                    elif inh == "Monoallelic":
                        self.monoallelic_homozygous_parents_filter(
                            v, variants[v]["child"], mum_gt, dad_gt, mum_aff, dad_aff
                        )
                    elif inh == "Mosaic":
                        self.mosaic_homozygous_parents_filter(v, variants[v]["child"], mum_gt, dad_gt, mum_aff, dad_aff)
                    elif inh == "Imprinted":
                        self.imprinted_homozygous_parents_filter(
                            v, variants[v]["child"], mum_gt, dad_gt, mum_aff, dad_aff
                        )
                    else:
                        logging.info(v + " unknown gene mode " + inh)
            else:
                logging.info(v + " fails inheritance filters - child must be " "hom or het")

    def autosomal_no_parents(self):
        """
        Screens variants in autosomes where there are both parents and adds
        candidates to candidate_variants
        """
        variants = self.variants_per_gene[self.hgncid]

        for v in variants.keys():
            if not variants[v]["child"].is_snv():
                continue
            for inh in self.gene["mode"]:
                if inh == "Biallelic":
                    self.biallelic_no_parents_filter(v, variants[v]["child"])
                elif inh == "Monoallelic":
                    self.monoallelic_no_parents_filter(v, variants[v]["child"])
                elif inh == "Mosaic":
                    self.mosaic_no_parents_filter(v, variants[v]["child"])
                elif inh == "Imprinted":
                    self.imprinted_no_parents_filter(v, variants[v]["child"])
                else:
                    logging.info(v + " unknown gene mode " + inh)

    def autosomal_single_parent(self):
        """
        Screens variants in autosomes where there is one parent and adds
        candidates to candidate_variants
        """
        pass

    def biallelic_heterozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        """
        Heterozygous variant in biallelic gene
        """
        self.inheritance_report.populate_inheritance_report(
            "autosomal", "biallelic", var.gt, mum_gt, dad_gt, mum_aff, dad_aff
        )

        vpass = "n"
        if mum_aff and dad_aff:
            if not (dad_gt == "1/1" and mum_gt == "1/1"):
                add_compound_het_to_candidates(varid, var, self.hgncid, "biallelic", self.candidate_variants)
                vpass = "y"
        elif mum_aff and not dad_aff:
            if not dad_gt == "1/1":
                add_compound_het_to_candidates(varid, var, self.hgncid, "biallelic", self.candidate_variants)
                vpass = "y"
        elif not mum_aff and dad_aff:
            if not mum_gt == "1/1":
                add_compound_het_to_candidates(varid, var, self.hgncid, "biallelic", self.candidate_variants)
                vpass = "y"
        else:
            if not dad_gt == "1/1" and not mum_gt == "1/1":
                add_compound_het_to_candidates(varid, var, self.hgncid, "biallelic", self.candidate_variants)
                vpass = "y"

        if vpass == "n":
            logging.info(varid + " failed inheritance filter for heterozygous " "variant in biallelic gene")

    def monoallelic_heterozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        """
        Heterozygous variant in monoallelic gene
        """
        self.inheritance_report.populate_inheritance_report(
            "autosomal", "monoallelic", var.gt, mum_gt, dad_gt, mum_aff, dad_aff
        )
        vpass = "n"

        if mum_aff and dad_aff:
            if not (dad_gt == "1/1" and mum_gt == "1/1"):
                add_single_var_to_candidates(varid, var, self.hgncid, "monoallelic", self.candidate_variants)
                vpass = "y"
        elif mum_aff and not dad_aff:
            if dad_gt == "0/0":
                add_single_var_to_candidates(varid, var, self.hgncid, "monoallelic", self.candidate_variants)
                vpass = "y"
        elif not mum_aff and dad_aff:
            if mum_gt == "0/0":
                add_single_var_to_candidates(varid, var, self.hgncid, "monoallelic", self.candidate_variants)
                vpass = "y"
        else:
            if mum_gt == "0/0" and dad_gt == "0/0":
                add_single_var_to_candidates(varid, var, self.hgncid, "monoallelic", self.candidate_variants)
                vpass = "y"

        if vpass == "n":
            logging.info(varid + " failed inheritance filter for heterozygous " "variant in monoallelic gene")

    def mosaic_heterozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        """
        Heterozygous variant in mosaic gene
        """
        self.inheritance_report.populate_inheritance_report(
            "autosomal", "mosaic", var.gt, mum_gt, dad_gt, mum_aff, dad_aff
        )

        vpass = "n"

        if mum_aff and dad_aff:
            if not (dad_gt == "1/1" and mum_gt == "1/1"):
                add_single_var_to_candidates(varid, var, self.hgncid, "mosaic", self.candidate_variants)
                vpass = "y"
        elif mum_aff and not dad_aff:
            if dad_gt == "0/0":
                add_single_var_to_candidates(varid, var, self.hgncid, "mosaic", self.candidate_variants)
                vpass = "y"
        elif not mum_aff and dad_aff:
            if mum_gt == "0/0":
                add_single_var_to_candidates(varid, var, self.hgncid, "mosaic", self.candidate_variants)
                vpass = "y"
        else:
            if mum_gt == "0/0" and dad_gt == "0/0":
                add_single_var_to_candidates(varid, var, self.hgncid, "mosaic", self.candidate_variants)
                vpass = "y"

        if vpass == "n":
            logging.info(varid + " failed inheritance filter for heterozygous " "variant in mosaic gene")

    def imprinted_heterozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        """
        Heterozygous variant in imprinted gene
        """
        self.inheritance_report.populate_inheritance_report(
            "autosomal", "imprinted", var.gt, mum_gt, dad_gt, mum_aff, dad_aff
        )
        vpass = "n"
        if mum_aff and dad_aff:
            if not (dad_gt == "1/1" and mum_gt == "1/1"):
                add_single_var_to_candidates(varid, var, self.hgncid, "imprinted", self.candidate_variants)
                vpass = "y"
        elif mum_aff and not dad_aff:
            if not dad_gt == "1/1":
                add_single_var_to_candidates(varid, var, self.hgncid, "imprinted", self.candidate_variants)
                vpass = "y"
        elif not mum_aff and dad_aff:
            if not mum_gt == "1/1":
                add_single_var_to_candidates(varid, var, self.hgncid, "imprinted", self.candidate_variants)
                vpass = "y"
        else:
            if not dad_gt == "1/1" and not mum_gt == "1/1":
                add_single_var_to_candidates(varid, var, self.hgncid, "imprinted", self.candidate_variants)
                vpass = "y"

        if vpass == "n":
            logging.info(varid + " failed inheritance filter for heterozygous " "variant in imprinted gene")

    def biallelic_homozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        """
        Homozygous variant in biallelic gene
        """
        self.inheritance_report.populate_inheritance_report(
            "autosomal", "biallelic", var.gt, mum_gt, dad_gt, mum_aff, dad_aff
        )

        vpass = "n"

        if mum_aff and dad_aff:
            if mum_gt != "0/0" and dad_gt != "0/0":
                add_single_var_to_candidates(varid, var, self.hgncid, "biallelic", self.candidate_variants)
                vpass = "y"
        elif mum_aff and not dad_aff:
            if dad_gt == "0/1" and mum_gt != "0/0":
                add_single_var_to_candidates(varid, var, self.hgncid, "biallelic", self.candidate_variants)
                vpass = "y"
        elif not mum_aff and dad_aff:
            if mum_gt == "0/1" and dad_gt != "0/0":
                add_single_var_to_candidates(varid, var, self.hgncid, "biallelic", self.candidate_variants)
                vpass = "y"
        else:
            if mum_gt == "0/1" and dad_gt == "0/1":
                add_single_var_to_candidates(varid, var, self.hgncid, "biallelic", self.candidate_variants)
                vpass = "y"
            elif mum_gt == "0/1" and dad_gt == "0/0":
                # allow homs apparently inherited from one parent as these may
                # have one hom and one cnv
                add_compound_het_to_candidates(varid, var, self.hgncid, "biallelic", self.candidate_variants)
                vpass = "y"
            elif mum_gt == "0/0" and dad_gt == "0/1":
                add_compound_het_to_candidates(varid, var, self.hgncid, "biallelic", self.candidate_variants)
                vpass = "y"

        if vpass == "n":
            logging.info(varid + " failed inheritance filter for homozygous " "variant in biallelic gene")

    def monoallelic_homozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        """
        Homozygous variant in monoallelic gene
        """
        # todo will need modification when CNVs (and UPDs) added
        self.inheritance_report.populate_inheritance_report(
            "autosomal", "monoallelic", var.gt, mum_gt, dad_gt, mum_aff, dad_aff
        )

        vpass = "n"

        if mum_aff and dad_aff:
            if mum_gt != "0/0" and dad_gt != "0/0":
                add_single_var_to_candidates(varid, var, self.hgncid, "monoallelic", self.candidate_variants)
                vpass = "y"

        if vpass == "n":
            logging.info(varid + " failed inheritance filter for homozygous " "variant in monoallelic gene")

    def mosaic_homozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        """
        Homozygous variant in mosaic gene
        """
        pass

    def imprinted_homozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        """
        Homozygous variant in imprinted gene
        """
        self.inheritance_report.populate_inheritance_report(
            "autosomal", "imprinted", var.gt, mum_gt, dad_gt, mum_aff, dad_aff
        )
        # all fail
        logging.info(varid + " failed inheritance filter for homozygous " "variant in imprinted gene")

    def biallelic_no_parents_filter(self, varid, var):
        """
        Variant in biallelic gene in a singleton
        """
        if var.genotype == "1":
            add_compound_het_to_candidates(varid, var, self.hgncid, "biallelic", self.candidate_variants)
        elif var.genotype == "2":
            add_single_var_to_candidates(varid, var, self.hgncid, "biallelic", self.candidate_variants)
        else:
            logging.info(varid + " failed inheritance filter for biallelic " "variant, invalid genotype")

    def monoallelic_no_parents_filter(self, varid, var):
        """
        Variant in monoalleic gene in a singleton
        """
        if var.genotype == "1" or var.genotype == "2":
            add_single_var_to_candidates(varid, var, self.hgncid, "monoallelic", self.candidate_variants)
        else:
            logging.info(varid + " failed inheritance filter for monoallelic " "variant, invalid genotype")

    def mosaic_no_parents_filter(self, varid, var):
        """
        Variant in mosaic gene in a singleton
        """
        if var.genotype == "1" or var.genotype == "2":
            add_single_var_to_candidates(varid, var, self.hgncid, "mosaic", self.candidate_variants)
        else:
            logging.info(varid + " failed inheritance filter for mosaic " "variant, invalid genotype")

    def imprinted_no_parents_filter(self, varid, var):
        """
        Variant in imprinted gene in a singleton
        """
        # todo will need modification when CNVs (and UPDs) added
        if var.genotype == "1":
            add_single_var_to_candidates(varid, var, self.hgncid, "imprinted", self.candidate_variants)
        elif var.genotype == "2":
            logging.info(varid + " failed inheritance filter for homozygous " "variant in imprinted gene")
            # todo add flag for CNV here or modify inheritance to include CNV data
        else:
            logging.info(varid + " failed inheritance filter for imprinted " "variant, invalid genotype")
