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

from utils.params import LOF_CONSEQUENCES, MODERATE_HIGH_IMPACT_CONSEQUENCES, REVEL_THRESHOLD_COMPOUNDHET_SINGLETON


class PostInheritanceFiltering(object):
    """
    Post-inheritance filters
    """

    def __init__(self, candidate_variants, family):
        self.candidate_variants = candidate_variants
        self.family = family

    def postinheritance_filter(self):
        """
        Post-inheritance filtering - MAF and REVEL
        """
        self.maf_filter()
        self.allele_count_filter()
        self.clean_spliceAI()
        return self.candidate_variants

    def maf_filter(self):
        """
        Filter non-Biallelic vairants with more stringent MAF thresholds
        """
        for v in list(self.candidate_variants["single_variants"].keys()):

            # We do not apply post-inheritance MAF filter to known pathogenic variants
            if (
                set(["Pathogenic", "Likely_pathogenic"])
                & set(self.candidate_variants["single_variants"][v]["variant"].ClinVar_CLNSIG.split("/"))
                != set()
            ):
                continue

            if "biallelic" not in self.candidate_variants["single_variants"][v]["mode"]:
                ddd_af = self.candidate_variants["single_variants"][v]["variant"].ddd_af
                max_af = self.candidate_variants["single_variants"][v]["variant"].max_af
                if ddd_af == ".":
                    ddd_af = "0"
                if max_af == ".":
                    max_af = "0"
                maximum_af = max(float(ddd_af), float(max_af))

                if self.family.has_both_parents() and maximum_af >= 0.0005:
                    del self.candidate_variants["single_variants"][v]
                    logging.info(
                        v + " failed post-inhertance MAF filter for family " "with parents, max AF = " + str(maximum_af)
                    )
                elif not self.family.has_both_parents() and maximum_af >= 0.0001:
                    del self.candidate_variants["single_variants"][v]
                    logging.info(
                        v + " failed post-inhertance MAF filter for family "
                        "without parents, max AF = " + str(maximum_af)
                    )

    def allele_count_filter(self):
        """
        Filter on AC_het and AC_hemi
        """
        for v in list(self.candidate_variants["single_variants"].keys()):

            # Not applicable to CNVs
            if not self.candidate_variants["single_variants"][v]["variant"].is_snv():
                continue

            if "biallelic" in self.candidate_variants["single_variants"][v]["mode"]:
                continue

            # We do not apply AC filters to known pathogenic variants
            if (
                set(["Pathogenic", "Likely_pathogenic"])
                & set(self.candidate_variants["single_variants"][v]["variant"].ClinVar_CLNSIG.split("/"))
                != set()
            ):
                continue

            if "monoallelic" in self.candidate_variants["single_variants"][v]["mode"]:
                if int(self.candidate_variants["single_variants"][v]["variant"].AC_het) > 12:
                    logging.info(
                        v + " failed post-inhertance AC_het filter for "
                        "monoallelic genes " + self.candidate_variants["single_variants"][v]["variant"].AC_het
                    )
                    del self.candidate_variants["single_variants"][v]
                    continue

            if ("hemizygous" in self.candidate_variants["single_variants"][v]["mode"]) and (
                self.candidate_variants["single_variants"][v]["variant"].sex == "XY"
            ):
                if int(self.candidate_variants["single_variants"][v]["variant"].AC_hemi) > 2:
                    logging.info(
                        v + " failed post-inhertance AC_hemi filter for "
                        "monoallelic genes " + self.candidate_variants["single_variants"][v]["variant"].AC_hemi
                    )
                    del self.candidate_variants["single_variants"][v]
                    continue

            if "X-linked dominant" in self.candidate_variants["single_variants"][v]["mode"]:
                AC_total = int(self.candidate_variants["single_variants"][v]["variant"].AC_het) + int(
                    self.candidate_variants["single_variants"][v]["variant"].AC_hemi
                )
                if AC_total > 12:
                    logging.info(
                        v + " failed post-inhertance AC_hemi filter for "
                        "X linked dominant genes "
                        + self.candidate_variants["single_variants"][v]["variant"].AC_het
                        + " + "
                        + self.candidate_variants["single_variants"][v]["variant"].AC_hemi
                    )
                    del self.candidate_variants["single_variants"][v]
                    continue

    def clean_spliceAI(self):
        """
        We want to keep high spliceAI variants only if they are :
        - single variant DNM
        - part of a compound het with the paired variant being LoF

        Here is the rationale for the filtering. The only reason why we would keep a low impact variant (e.g. synonymous)
        or a missense variant with low REVEL score, is because it has a high spliceAI score.

        #TODO : we should normally doublecheck the spliceAI score, however we will implement this in the next version of the tool

        """

        # Single variant filter
        for varid in list(self.candidate_variants["single_variants"].keys()):
            variant = self.candidate_variants["single_variants"][varid]["variant"]

            # We keep single variants DNM in any case
            if variant.is_snv() and variant.is_dnm():
                continue

            # If a missense variant passed the pre-inheritance filters because of high spliceAI score but is not part of a compound het,
            # we now remove it
            if variant.consequence == "missense_variant":
                if variant.revel == "." or float(variant.revel) < 0.4:
                    logging.info(varid + " failed post-inheritance spliceAI filter")
                    del self.candidate_variants["single_variants"][varid]

            # If a variant passed the pre-inheritance filters because of high spliceAI score but is not part of a compound het,
            # we now remove it
            if variant.consequence not in MODERATE_HIGH_IMPACT_CONSEQUENCES:
                logging.info(varid + " failed post-inheritance spliceAI filter")
                del self.candidate_variants["single_variants"][varid]

        # Compound het filter
        for compound_het_id in list(self.candidate_variants["compound_hets"].keys()):

            compound_het = self.candidate_variants["compound_hets"][compound_het_id]
            assert len(compound_het) == 2
            var1 = compound_het[list(compound_het.keys())[0]]["variant"]
            var2 = compound_het[list(compound_het.keys())[1]]["variant"]

            # Compound het with both variants missene or equivalent with a high REVEL score or known pathogenic in ClinVar are kept
            if (var1.revel != "." and float(var1.revel) > REVEL_THRESHOLD_COMPOUNDHET_SINGLETON) or (
                set(["Pathogenic", "Likely_pathogenic"]) & set(var1.ClinVar_CLNSIG.split("/")) != set()
            ):

                if (var2.revel != "." and float(var2.revel) > REVEL_THRESHOLD_COMPOUNDHET_SINGLETON) or (
                    set(["Pathogenic", "Likely_pathogenic"]) & set(var2.ClinVar_CLNSIG.split("/")) != set()
                ):
                    continue

            # If one variant has a high spliceAI score but the other one is not a LoF we remove the compound het
            if (var1.consequence not in MODERATE_HIGH_IMPACT_CONSEQUENCES) or (
                (var1.consequence == "missense_variant") and (var1.revel == "." or float(var1.revel) < 0.4)
            ):
                if var2.consequence not in LOF_CONSEQUENCES:
                    logging.info(compound_het_id + " failed post-inheritance spliceAI filter for compound hets")
                    del self.candidate_variants["compound_hets"][compound_het_id]
                    continue

            # If one variant has a high spliceAI score but the other one is not a LoF we remove the compound het
            if (var2.consequence not in MODERATE_HIGH_IMPACT_CONSEQUENCES) or (
                (var2.consequence == "missense_variant") and (var2.revel == "." or float(var2.revel) < 0.4)
            ):

                if var1.consequence not in LOF_CONSEQUENCES:
                    logging.info(compound_het_id + " failed post-inheritance spliceAI filter for compound hets")
                    del self.candidate_variants["compound_hets"][compound_het_id]
                    continue
