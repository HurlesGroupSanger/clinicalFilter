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

from utils.utils import add_single_var_to_candidates
from utils.utils import convert_genotype_to_gt


class AllosomalFilter(object):
    """
    Inheritance filters for SNVs/indels on allosomes
    """

    def __init__(self, Inheritancefilter, hgncid):
        self.family = Inheritancefilter.family
        self.parents = Inheritancefilter.parents
        self.variants_per_gene = Inheritancefilter.variants_per_gene
        self.candidate_variants = Inheritancefilter.candidate_variants
        self.inheritance_report = Inheritancefilter.inhreport
        self.gene = Inheritancefilter.genes[hgncid]
        self.hgncid = hgncid
        self.proband_X_count = self.family.proband.X_count

    def allosomal_filter(self):
        if self.parents == "both":
            self.allosomal_both_parents()
        elif self.parents == "none":
            self.allosomal_no_parents()
        else:
            self.allosomal_single_parent()

    def allosomal_both_parents(self):
        """
        screens variants in X/Y where there are both parents and adds candidates
        to candidate_variants
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

            genotype = self.get_variant_genotype(variants[v]["child"], v)
            # if dad gt = 0/1 should go to a different (mosaic) pipeline -
            # fail for now
            # TODO : not sure about this has it does not even do any action, just logging
            if dad_gt == "0/1" and variants[v]["child"].chrom == "X":
                logging.info(v + " failed due to 0/1 paternal genotype in X: ")
            # if genotype is none then variant fails
            if genotype is None:
                continue
            elif genotype == "homozygous":

                # Variants found on both copies of X chromosome in female

                for inh in self.gene["mode"]:
                    # No such thing as homozygous variant in hemizygous gene
                    if inh == "Hemizygous":
                        self.gn_hemizygous_gt_homozygous_parents_filter(
                            v, variants[v]["child"], mum_gt, dad_gt, mum_aff, dad_aff
                        )
                    # TODO : not sure why no variants are kept here
                    # Is it considered to unlikely to have an inherited hom ?
                    elif inh == "X-linked dominant":
                        self.gn_X_linked_dominant_gt_homozygous_parents_filter(
                            v, variants[v]["child"], mum_gt, dad_gt, mum_aff, dad_aff
                        )
                    # TODO : not sure why no variants are kept here
                    # Is it considered to unlikely to have an inherited hom ?
                    elif inh == "X-linked over-dominance":
                        self.gn_X_linked_over_dominant_gt_homozygous_parents_filter(
                            v, variants[v]["child"], mum_gt, dad_gt, mum_aff, dad_aff
                        )
                    # No such thing as homozygous variant in Y gene
                    elif inh == "monoallelic_Y_hem":
                        logging.info(v + "fails inheritance filters homozygous GT in " + inh)
                    else:
                        logging.info(v + " unknown gene mode " + inh)
            elif genotype == "hemizygous":

                # Variants found on X or Y chromosome in male

                for inh in self.gene["mode"]:
                    # Hemyzygous variant found in Hemyzygous gene
                    # Keep variant unless the mother has two copies but is unaffected
                    if inh == "Hemizygous":
                        self.gn_hemizygous_gt_hemizygous_parents_filter(
                            v, variants[v]["child"], mum_gt, dad_gt, mum_aff, dad_aff
                        )
                    # Hemyzygous variant found in X-linked dominant gene
                    # Keep variant unless the mother has two copies but is unaffected
                    elif inh == "X-linked dominant":
                        self.gn_X_linked_dominant_gt_hemizygous_parents_filter(
                            v, variants[v]["child"], mum_gt, dad_gt, mum_aff, dad_aff
                        )
                    # Hemyzygous variant found in X-linked over-dominance gene
                    # Discard variant. TODO : Why ?
                    elif inh == "X-linked over-dominance":
                        self.gn_X_linked_over_dominant_gt_hemizygous_parents_filter(
                            v, variants[v]["child"], mum_gt, dad_gt, mum_aff, dad_aff
                        )
                    # Hemyzygous variant found in monoallelic_Y_hem gene
                    # Keep variant unless the father has the variant but is unaffected
                    elif inh == "monoallelic_Y_hem":
                        self.gn_mono_y_hem_gt_hemizygous_parents_filter(
                            v, variants[v]["child"], mum_gt, dad_gt, mum_aff, dad_aff
                        )
                    else:
                        logging.info(v + " unknown gene mode " + inh)
            elif genotype == "heterozygous":

                # Variants found on a single copy of X chromosome in female

                for inh in self.gene["mode"]:

                    # Heterozygous variant found in Hemizygous gene
                    # Discard variant if found in unaffected mother or father
                    if inh == "Hemizygous":
                        self.gn_hemizygous_gt_heterozygous_parents_filter(
                            v, variants[v]["child"], mum_gt, dad_gt, mum_aff, dad_aff
                        )

                    # Heterozygous variant found in X-linked dominant gene
                    # Discard variant if found in unaffected mother or father
                    elif inh == "X-linked dominant":
                        self.gn_X_linked_dominant_gt_heterozygous_parents_filter(
                            v, variants[v]["child"], mum_gt, dad_gt, mum_aff, dad_aff
                        )

                    # Heterozygous variant found in X-linked over-dominance gene
                    # Keep variant if the variant is denovo, or :
                    # - if the father is unaffected, the mother affected and has a single copy of the variant
                    # - if the father is unaffected, the mother unaffected and has no or two copies of the variant
                    elif inh == "X-linked over-dominance":
                        self.gn_X_linked_over_dominant_gt_heterozygous_parents_filter(
                            v, variants[v]["child"], mum_gt, dad_gt, mum_aff, dad_aff
                        )
                    # No such thing as heterozygous variant in Y gene
                    elif inh == "monoallelic_Y_hem":
                        logging.info(v + "fails inheritance filters heterozygous GT in " + inh)
                    else:
                        logging.info(v + " unknown gene mode " + inh)
            else:
                logging.info(v + " unknown genotype " + genotype)

    def allosomal_no_parents(self):
        """
        Screens variants in X/Y where there are no parents and adds candidates
        to candidate_variants
        """
        # All allosomal variants on X with no parents pass, Y fail if 0/1
        variants = self.variants_per_gene[self.hgncid]
        for v in variants.keys():
            if not variants[v]["child"].is_snv():
                continue
            genotype = self.get_variant_genotype(variants[v]["child"], v)
            if genotype is None:
                continue
            if variants[v]["child"].chrom == "X":
                for inh in self.gene["mode"]:
                    add_single_var_to_candidates(
                        v,
                        variants[v]["child"],
                        self.hgncid,
                        inh.lower(),
                        self.candidate_variants,
                    )
            elif variants[v]["child"].chrom == "Y":
                for inh in self.gene["mode"]:
                    if variants[v]["child"].gt == "1/1":
                        add_single_var_to_candidates(
                            v,
                            variants[v]["child"],
                            self.hgncid,
                            inh.lower(),
                            self.candidate_variants,
                        )
                    else:
                        logging.info(v + "fails inheritance filters non-hemizygous GT " "in " + inh)
            else:
                logging.info(
                    v + "fails inheritance allosomal inheritance filters, " "chromosome = " + variants[v]["child"].chrom
                )

    def allosomal_single_parent(self):
        """
        Screens variants in X/Y where there is one parent and adds candidates
        to candidate_variants
        """
        pass

    def gn_hemizygous_gt_homozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        """
        hemizygous gene, homozygous proband all fail
        """
        self.inheritance_report.populate_inheritance_report(
            "allosomal", "hemizygous", "homozygous", mum_gt, dad_gt, mum_aff, dad_aff
        )
        logging.info(varid + " failed inheritance filter for homozygous " "variant in hemizygous gene")

    def gn_X_linked_dominant_gt_homozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        """
        X-linked dominant gene, homozygous proband all fail
        """
        self.inheritance_report.populate_inheritance_report(
            "allosomal",
            "X-linked_dominant",
            "homozygous",
            mum_gt,
            dad_gt,
            mum_aff,
            dad_aff,
        )
        logging.info(varid + " failed inheritance filter for homozygous " "variant in X-linked dominant gene")

    def gn_X_linked_over_dominant_gt_homozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        """
        X-linked over dominant gene, homozygous proband all fail
        """
        self.inheritance_report.populate_inheritance_report(
            "allosomal",
            "X-linked_over_dominance",
            "homozygous",
            mum_gt,
            dad_gt,
            mum_aff,
            dad_aff,
        )
        logging.info(varid + " failed inheritance filter for homozygous " "variant in X-linked over dominance gene")

    def gn_hemizygous_gt_hemizygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        """
        hemizygous gene, hemizygous proband pass unless mum is 1/1 and
        unaffected
        """
        self.inheritance_report.populate_inheritance_report(
            "allosomal", "hemizygous", "hemizygous", mum_gt, dad_gt, mum_aff, dad_aff
        )
        # If the mother has two copies of the variant but is unaffected, then discard the variant
        if not mum_aff and mum_gt == "1/1":
            logging.info(varid + " failed inheritance filter for hemizygous " "variant in hemizygous gene")
        # Otherwise consider it a plausible candidate
        else:
            add_single_var_to_candidates(varid, var, self.hgncid, "hemizygous", self.candidate_variants)

    def gn_X_linked_dominant_gt_hemizygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        """
        X linked dominant gene and hemizygous proband pass unless mum is 1/1
        and unaffected
        """
        self.inheritance_report.populate_inheritance_report(
            "allosomal",
            "X-linked_dominant",
            "hemizygous",
            mum_gt,
            dad_gt,
            mum_aff,
            dad_aff,
        )

        # If the mother has two copies of the variant but is unaffected, then discard the variant
        if not mum_aff and mum_gt == "1/1":
            logging.info(varid + " failed inheritance filter for hemizygous " "variant in X-linked dominant gene")
        # Otherwise consider it a plausible candidate
        else:
            add_single_var_to_candidates(varid, var, self.hgncid, "X-linked dominant", self.candidate_variants)

    def gn_X_linked_over_dominant_gt_hemizygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        """
        X-linked over dominant gene, hemizygous proband all fail
        """
        self.inheritance_report.populate_inheritance_report(
            "allosomal",
            "X-linked_over_dominance",
            "hemizygous",
            mum_gt,
            dad_gt,
            mum_aff,
            dad_aff,
        )
        logging.info(varid + " failed inheritance filter for hemizygous " "variant in X-linked over dominance gene")

    def gn_hemizygous_gt_heterozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        """
        Hemizygous gene, heterozygous proband. Fail if dad is 1/1 unaffected or
        mum 0/1 or 1/1 and unaffected
        """
        self.inheritance_report.populate_inheritance_report(
            "allosomal", "hemizygous", "heterozygous", mum_gt, dad_gt, mum_aff, dad_aff
        )

        # If the mother has at least one copy of the variant but is unaffected, discard variant
        if not mum_aff and mum_gt != "0/0":
            logging.info(varid + " failed inheritance filter for heterozygous " "variant in hemizygous gene")

        # If the father has the variant but is unaffected, discard variant
        elif not dad_aff and dad_gt != "0/0":
            logging.info(varid + " failed inheritance filter for heterozygous " "variant in hemizygous gene")
        # Otherwise consider it a plausible candidate
        else:
            add_single_var_to_candidates(varid, var, self.hgncid, "hemizygous", self.candidate_variants)

    def gn_X_linked_dominant_gt_heterozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        """
        X-linked dominant gene, heterozygous proband. Fail if dad is 1/1 unaffected or
        mum 0/1 or 1/1 and unaffected
        """
        self.inheritance_report.populate_inheritance_report(
            "allosomal",
            "X-linked_dominant",
            "heterozygous",
            mum_gt,
            dad_gt,
            mum_aff,
            dad_aff,
        )

        # If the mother has at least one copy of the variant but is unaffected, discard variant
        if not mum_aff and mum_gt != "0/0":
            logging.info(varid + " failed inheritance filter for heterozygous " "variant in X-linked dominant gene")
        # If the father has the variant but is unaffected, discard variant
        elif not dad_aff and dad_gt != "0/0":
            logging.info(varid + " failed inheritance filter for heterozygous " "variant in X-linked dominant gene")
        # Otherwise consider it a plausible candidate
        else:
            add_single_var_to_candidates(varid, var, self.hgncid, "X-linked dominant", self.candidate_variants)

    def gn_X_linked_over_dominant_gt_heterozygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        """
        X-linked over-dominant gene, heterozygous proband. If dad aff pass if
        DNM if dad unaff pass if DNM or mum_aff and 1/1 or mum unaff and not 0/1

        Overdominance is a genetic phenomenon that occurs when a heterozygote's phenotype
        is more extreme or better adapted than either of its homozygous parents
        """
        self.inheritance_report.populate_inheritance_report(
            "allosomal",
            "X-linked_over_dominance",
            "heterozygous",
            mum_gt,
            dad_gt,
            mum_aff,
            dad_aff,
        )

        # Keep variant if father affected and variant is denovo
        if dad_aff:
            if mum_gt == "0/0" and dad_gt == "0/0":
                add_single_var_to_candidates(
                    varid,
                    var,
                    self.hgncid,
                    "X-linked over-dominance",
                    self.candidate_variants,
                )
            else:
                logging.info(
                    varid + " failed inheritance filter for heterozygous " "variant in X-linked over-dominance gene"
                )
        else:

            # Keep variant if father unaffected and variant is denovo
            if mum_gt == "0/0" and dad_gt == "0/0":
                add_single_var_to_candidates(
                    varid,
                    var,
                    self.hgncid,
                    "X-linked over-dominance",
                    self.candidate_variants,
                )

            # Keep variant if father unaffected, mother affected and has a single copy of the variant
            elif mum_aff and mum_gt == "0/1":
                add_single_var_to_candidates(
                    varid,
                    var,
                    self.hgncid,
                    "X-linked over-dominance",
                    self.candidate_variants,
                )
            # Keep variant if father unaffected, mother unaffected and has no or two copies of the variant
            elif not mum_aff and not mum_gt == "0/1":
                add_single_var_to_candidates(
                    varid,
                    var,
                    self.hgncid,
                    "X-linked over-dominance",
                    self.candidate_variants,
                )
            else:
                logging.info(
                    varid + " failed inheritance filter for heterozygous " "variant in X-linked over-dominance gene"
                )

    def gn_mono_y_hem_gt_hemizygous_parents_filter(self, varid, var, mum_gt, dad_gt, mum_aff, dad_aff):
        """
        Fail unless dad affected and 1/1
        """
        self.inheritance_report.populate_inheritance_report(
            "allosomal",
            "monoallelic_Y_hemizygous",
            "hemizygous",
            mum_gt,
            dad_gt,
            mum_aff,
            dad_aff,
        )

        # If the father has the variant but is unaffected, then discard the variant
        # as the father should be affected
        if not dad_aff and dad_gt == "1/1":
            logging.info(
                varid + " failed inheritance filter for hemizygous " "variant in monoallelic_Y_hemizygous gene"
            )

        # Otherwise consider it a plausible candidate
        else:
            add_single_var_to_candidates(
                varid,
                var,
                self.hgncid,
                "monoallelic_Y_hemizygous",
                self.candidate_variants,
            )

    def get_variant_genotype(self, variant, v):
        """
        Work out if a variant is homozygous, heterozygous or hemizygous
        """
        if variant.chrom == "X":
            if self.proband_X_count >= 2:
                if variant.gt == "0/1":
                    genotype = "heterozygous"
                elif variant.gt == "1/1":
                    genotype = "homozygous"
                else:
                    logging.info(v + " fails invalid genotype " + variant.gt)
                    return None
            elif self.proband_X_count == 1:
                if variant.gt == "1/1":
                    genotype = "hemizygous"
                elif variant.gt == "0/1":
                    adsplit = variant.ad.split(",")
                    vaf = int(adsplit[1]) / (int(adsplit[0]) + int(adsplit[1]))
                    if vaf > 0.8 or variant.dnm == True:
                        genotype = "hemizygous"
                    else:
                        logging.info(
                            v + " fails 0/1 variant with low VAF in proband " "with 1 X chromosome and not DNM"
                        )
                        return None
                else:
                    logging.info(v + " fails invalid genotype " + variant.gt)
                    return None
            else:
                logging.info(v + " fails invalid X chromsome count")
                return None
        elif variant.chrom == "Y":
            # TODO : how can we have heterozygous variant on Y ?
            if variant.gt == "0/1":
                genotype = "heterozygous"
            elif variant.gt == "1/1":
                genotype = "hemizygous"
            else:
                logging.info(v + " fails invalid genotype " + variant.gt)
                return None

        return genotype
