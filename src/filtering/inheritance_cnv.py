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
from utils.utils import add_compound_het_to_candidates


class CNVFiltering(object):
    """
    Inheritance filters for CNVs
    """

    def __init__(self, variants, family, genes, regions, trusted_variants, candidate_variants):
        self.variants = variants
        self.family = family
        self.genes = genes
        self.regions = regions
        self.trusted_variants = trusted_variants
        self.candidate_variants = candidate_variants
        self.parents = None
        if self.family.has_both_parents():
            self.parents = "both"
        elif self.family.has_no_parents():
            self.parents = "none"
        elif self.family.has_dad():
            self.parents = "dad_only"
        elif self.family.has_mum():
            self.parents = "mum_only"
        if self.family.has_mum():
            self.mum_aff = self.family.mum.affected
        if self.family.has_dad():
            self.dad_aff = self.family.dad.affected

    def cnv_filter(self):
        if self.genes:
            self.cnv_filter_genes()
        if self.regions:
            # todo filtering for regions
            pass
        if self.trusted_variants:
            # todo filtering for trusted variants
            pass

    def cnv_filter_genes(self):
        if self.parents == "both":
            self.cnv_filter_parents()
        elif self.parents == "none":
            self.cnv_filter_no_parents()
        else:
            self.cnv_filter_single_parent()

    def cnv_filter_parents(self):
        for v in self.variants["child"].keys():
            if not self.variants["child"][v].is_cnv():
                continue
            if self.variants["child"][v].cnv_filter == "Fail":
                logging.info(v + " CNV failed quality filters")
                continue
            if self.genes:
                modes = self.get_ddg2p_modes(v)
            # inheritance
            self.inhmatch = self.cnv_inheritance_filter(v)
            if not self.inhmatch:
                # possible compound het
                self.posscomphet = self.cnv_candidate_compound_het_filter(v, modes)
                if not self.posscomphet:
                    logging.info(v + " failed CNV filter, inheritance doesn't " "match and not possible compound het")
            else:
                # non-ddg2p filter
                self.passnonddg2p = self.cnv_non_ddg2p_filter(v)
                if not self.passnonddg2p:
                    # ddg2p filter
                    self.passddg2p = self.cnv_ddg2p_filter(v)
                    if not self.passddg2p:
                        self.posscomphet = self.cnv_candidate_compound_het_filter(v, modes)
                        if not self.posscomphet:
                            logging.info(
                                v + " failed CNV filter, inheritance doesn't " "match and not possible compound het"
                            )

    def cnv_filter_no_parents(self):
        # non-ddg2p filter
        for v in self.variants["child"].keys():
            if not self.variants["child"][v].is_cnv():
                continue
            if self.variants["child"][v].cnv_filter == "Fail":
                logging.info(v + " CNV failed quality filters")
                continue
            if self.genes:
                modes = self.get_ddg2p_modes(v)
            self.passnonddg2p = self.cnv_non_ddg2p_filter(v)
            if not self.passnonddg2p:
                # ddg2p filter
                self.passddg2p = self.cnv_ddg2p_filter(v)
                if not self.passddg2p:
                    self.posscomphet = self.cnv_candidate_compound_het_filter(v, modes)
                    if not self.posscomphet:
                        logging.info(v + " failed CNV filter and not possible compound " "het")

    def cnv_filter_single_parent(self):
        pass

    def cnv_inheritance_filter(self, varid):
        """
        CNV inheritance
        return True or False for pass or fail
        inheritance matches parental affected status if:
        paternal and father affected
        OR
        maternal and mum affected
        OR
        biparental plus CN = 0
        OR
        biparental and both/either parent(s) affected
        OR
        male (XY) proband X chromosome, maternal inh and mum unaffected
        no variants are added to candidates at this stage
        """
        if not self.variants["child"][varid].cnv_inh in [
            "paternal_inh",
            "maternal_inh",
            "biparental_inh",
        ]:
            return True
        if self.variants["child"][varid].cnv_inh == "paternal_inh" and self.dad_aff:
            return True
        elif self.variants["child"][varid].cnv_inh == "maternal_inh" and self.mum_aff:
            return True
        elif self.variants["child"][varid].cnv_inh == "biparental_inh" and self.variants["child"][varid].cn == "0":
            return True
        elif self.variants["child"][varid].cnv_inh == "biparental_inh" and (self.mum_aff or self.dad_aff):
            return True
        elif (
            self.family.proband.sex == "XY"
            and self.variants["child"][varid].chrom == "X"
            and not self.mum_aff
            and self.variants["child"][varid].cnv_inh == "maternal_inh"
        ):
            return True
        else:
            return False

    def cnv_candidate_compound_het_filter(self, varid, modes):
        """
        Identify CNVS which could be in compound hets
        """

        # TODO : why desired_cn is 1 or 3 ???
        desired_cn = ["1", "3"]

        # could the CNV be part of a compound het? If so, add to candidate
        # compound hets
        if not self.variants["child"][varid].cn in desired_cn:
            return False
        else:
            # is any gene covered by the CNV biallelic, or CN = 1 and male
            # hemizygous
            if "Biallelic" in modes.keys():
                for hgncid in modes["Biallelic"]:
                    self.variants["child"][varid].reportable_symbol.append(self.genes[hgncid]["symbol"])
                    self.variants["child"][varid].reportable_hgnc_id.append(hgncid)
                    add_compound_het_to_candidates(
                        varid,
                        self.variants["child"][varid],
                        hgncid,
                        "Biallelic",
                        self.candidate_variants,
                    )
                return True
            elif (
                "Hemizygous" in modes.keys()
                and self.variants["child"][varid].cn == "1"
                and self.family.proband.sex == "XY"
            ):
                for hgncid in modes["Hemizygous"]:
                    self.variants["child"][varid].reportable_symbol.append(self.genes[hgncid]["symbol"])
                    self.variants["child"][varid].reportable_hgnc_id.append(hgncid)
                    add_compound_het_to_candidates(
                        varid,
                        self.variants["child"][varid],
                        hgncid,
                        "Hemizygous",
                        self.candidate_variants,
                    )
                return True
            else:
                return False

    def cnv_non_ddg2p_filter(self, varid):
        """
        CNVs of >1M pass regardless of gene content
        """
        if int(self.variants["child"][varid].cnv_length) > 1000000:
            add_single_var_to_candidates(varid, self.variants["child"][varid], "-", "-", self.candidate_variants)
            return True
        else:
            return False

    def cnv_ddg2p_filter(self, varid):
        """
        CNV DDG2P filter
        """
        cnvpass = False
        # get all genes covered by the CNV and go through each one at a time to
        # see if any pass
        hgncids = self.variants["child"][varid].hgnc_id_all.split("|")
        for hid in hgncids:

            surrounding_dup = False

            # Extract HGNC number only (e.g., HGNC:114 will return 114) and check if in DDG2P genes
            hgncid = hid[5:]
            if not hgncid in self.genes.keys():
                continue

            # fail duplications completely surrounding surround monoallelic,
            # hemizygous and x-linked dominant genes with loss of function
            # mechanism
            if (
                self.variants["child"][varid].alt == "<DUP>"
                and self.variants["child"][varid].chrom == self.genes[hgncid]["chr"]
            ):
                dupmodes = set({"Monoallelic", "Hemizygous", "X-linked dominant"})

                if (
                    "Loss of function" in self.genes[hgncid]["mechanism"]
                    and len(set.intersection(self.genes[hgncid]["mode"], dupmodes)) > 0
                ):
                    if (int(self.variants["child"][varid].pos) < int(self.genes[hgncid]["start"])) and (
                        int(self.variants["child"][varid].cnv_end) > int(self.genes[hgncid]["end"])
                    ):
                        logging.info(
                            varid + " duplication completely surrounds "
                            "monoallelic, hemizygous or heterozygus "
                            "gene with LoF mechanism " + self.genes[hgncid]["symbol"] + " hgnc:" + hgncid
                        )
                        surrounding_dup = True

            # Biallelic gene pass if copy number (CN) = 0 and mechanism in
            # "Uncertain", "Loss of function", "Dominant negative"
            if int(self.variants["child"][varid].cn) == 0 and "Biallelic" in self.genes[hgncid]["mode"]:
                biallelicmechs = set({"Uncertain", "Loss of function", "Dominant negative"})
                if len(set.intersection(self.genes[hgncid]["mechanism"], biallelicmechs)) > 0:
                    cnvpass = True
                    self.variants["child"][varid].reportable_symbol.append(self.genes[hgncid]["symbol"])
                    self.variants["child"][varid].reportable_hgnc_id.append(hgncid)
                    add_single_var_to_candidates(
                        varid,
                        self.variants["child"][varid],
                        hgncid,
                        "biallelic",
                        self.candidate_variants,
                    )
                    return cnvpass

            # Monoallelic, X-linked dominant or Hemizygous in male pass if
            # CN=0, 1 or 3 and any mechanism
            cns_wanted = ["0", "1", "3"]
            if not surrounding_dup:
                if "Monoallelic" in self.genes[hgncid]["mode"] or "X-linked dominant" in self.genes[hgncid]["mode"]:
                    if self.variants["child"][varid].cn in cns_wanted:
                        self.variants["child"][varid].reportable_symbol.append(self.genes[hgncid]["symbol"])
                        self.variants["child"][varid].reportable_hgnc_id.append(hgncid)
                        cnvpass = True
                        add_single_var_to_candidates(
                            varid,
                            self.variants["child"][varid],
                            hgncid,
                            (",").join(self.genes[hgncid]["mode"]),
                            self.candidate_variants,
                        )
                        return cnvpass
                if "Hemizygous" in self.genes[hgncid]["mode"] and self.family.proband.sex == "XY":
                    cnvpass = True
                    self.variants["child"][varid].reportable_symbol.append(self.genes[hgncid]["symbol"])
                    self.variants["child"][varid].reportable_hgnc_id.append(hgncid)
                    add_single_var_to_candidates(
                        varid,
                        self.variants["child"][varid],
                        hgncid,
                        (",").join(self.genes[hgncid]["mode"]),
                        self.candidate_variants,
                    )
                    return cnvpass

                # Hemizygous in female pass if CN=3 and mechanism =
                # "Increased gene dosage"
                if (
                    "Hemizygous" in self.genes[hgncid]["mode"]
                    and self.family.proband.sex == "XX"
                    and "Increased gene dosage" in self.genes[hgncid]["mechanism"]
                    and self.variants["child"][varid].cn == "3"
                ):
                    cnvpass = True
                    self.variants["child"][varid].reportable_symbol.append(self.genes[hgncid]["symbol"])
                    self.variants["child"][varid].reportable_hgnc_id.append(hgncid)
                    add_single_var_to_candidates(
                        varid,
                        self.variants["child"][varid],
                        hgncid,
                        "Hemizygous",
                        self.candidate_variants,
                    )
                    return cnvpass

            # Pass intragenic DUP in monoallelic or X-linked dominant gene with
            # loss of function mechanism and any part of the gene is outside of
            # the CNV boundary
            if self.variants["child"][varid].alt == "<DUP>":
                if (
                    "Monoallelic" in self.genes[hgncid]["mode"] or "X-linked dominant" in self.genes[hgncid]["mode"]
                ) and "Loss of function" in self.genes[hgncid]["mechanism"]:
                    if int(self.variants["child"][varid].pos) > int(self.genes[hgncid]["start"]) or int(
                        self.variants["child"][varid].cnv_end
                    ) < int(self.genes[hgncid]["end"]):
                        cnvpass = True
                        self.variants["child"][varid].reportable_symbol.append(self.genes[hgncid]["symbol"])
                        self.variants["child"][varid].reportable_hgnc_id.append(hgncid)
                        add_single_var_to_candidates(
                            varid,
                            self.variants["child"][varid],
                            hgncid,
                            (",").join(self.genes[hgncid]["mode"]),
                            self.candidate_variants,
                        )
                        return cnvpass

        return cnvpass

    def get_ddg2p_modes(self, varid):
        modes = {}
        # when a gene list is given, find the modes of all
        hgncids = self.variants["child"][varid].hgnc_id_all.split("|")
        for hid in hgncids:
            hgncid = hid[5:]
            if hgncid in self.genes.keys():
                genemodes = self.genes[hgncid]["mode"]
                for m in genemodes:
                    if m not in modes.keys():
                        modes[m] = []
                    modes[m].append(hgncid)
        return modes
