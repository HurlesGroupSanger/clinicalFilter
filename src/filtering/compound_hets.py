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
from itertools import combinations
from utils.utils import common_elements


class CompoundHetScreen(object):
    """
    Class for screening variants and identifying compound hets
    """

    def __init__(self, candidate_variants, family):
        self.candidate_variants = candidate_variants
        self.family = family

    def screen_compound_hets(self):
        # sort out compound hets
        compound_het_passes = {}
        for gn in self.candidate_variants["compound_hets"].keys():
            # print(gn)
            if len(self.candidate_variants["compound_hets"][gn].keys()) < 2:
                for v in self.candidate_variants["compound_hets"][gn].keys():
                    logging.info(v + " failed compound het screen: <2 vars in hgnc " + gn)

            combs_to_screen = list(combinations(self.candidate_variants["compound_hets"][gn].keys(), 2))
            for pair in combs_to_screen:
                var1 = self.candidate_variants["compound_hets"][gn][pair[0]]["variant"]
                var2 = self.candidate_variants["compound_hets"][gn][pair[1]]["variant"]

                if var1 != var2:
                    if self.is_compound_pair(pair[0], var1, pair[1], var2):
                        if not gn in compound_het_passes.keys():
                            compound_het_passes[gn] = {}
                        if not pair[0] in compound_het_passes[gn].keys():
                            compound_het_passes[gn][pair[0]] = {}
                        if not pair[1] in compound_het_passes[gn].keys():
                            compound_het_passes[gn][pair[1]] = {}

                        compound_het_passes[gn][pair[0]]["variant"] = var1
                        compound_het_passes[gn][pair[0]]["mode"] = self.candidate_variants["compound_hets"][gn][
                            pair[0]
                        ]["mode"]
                        compound_het_passes[gn][pair[1]]["variant"] = var2
                        compound_het_passes[gn][pair[1]]["mode"] = self.candidate_variants["compound_hets"][gn][
                            pair[1]
                        ]["mode"]

        self.candidate_variants["compound_hets"] = compound_het_passes

    def is_compound_pair(self, varid1, var1, varid2, var2):
        """
        Test to see if a pair of variants could be a compound het
        """

        if self.family.has_no_parents():
            # If there are no parents and the variants are not both missense
            # or in frame del/ins then pass
            missense_equiv = ["missense_variant", "inframe_deletion", "inframe_insertion"]
            var1cqs = var1.consequence.split("&")
            var2cqs = var2.consequence.split("&")
            var1_missense_equiv = common_elements(missense_equiv, var1cqs)
            var2_missense_equiv = common_elements(missense_equiv, var2cqs)
            if len(var1_missense_equiv) > 0 and len(var2_missense_equiv) > 0:
                # if var1.consequence.find(
                #         "missense_variant") and var2.consequence.find(
                #     "missense_variant"):
                logging.info(
                    varid1 + " " + varid2 + " failed compound het screen, no parents and " "both missense or equivalent"
                )
                return False
            elif (var1.pid != "." and var2.pid != ".") and (var1.pid == var2.pid):
                logging.info(varid1 + " " + varid2 + " failed compound het " "screen, variants in cis")
                return False

            else:
                return True
        elif self.family.has_both_parents():
            if (
                var1.chrom == "X"
                and not self.family.dad.get_affected_status()
                and (var1.get_dad_genotype() == "0" or var2.get_dad_genotype() == "0")
            ):
                logging.info(
                    varid1 + " " + varid2 + " failed compound het screen, X chrom and dad "
                    "unaffected and hom ref for 1 variant"
                )
                return False
            elif var1.triogenotype in ["201", "210"] or var1.triogenotype in ["201", "210"]:
                # triogenotype of 201 or 210 only passes is one variant is a
                # deletion with copy number 1
                if var1.is_cnv() and not var2.is_cnv():
                    varcnv = var1
                    varsnv = var2
                elif var2.is_cnv() and not var1.is_cnv():
                    varcnv = var2
                    varsnv = var1
                else:
                    logging.info(
                        varid1 + " " + varid2 + " failed compound het screen, homozygous on "
                        "one alllle and other allele not CNV"
                    )
                if varcnv.cn == "1":
                    if varsnv.triogenotype == "201" and varcnv.triogenotype == "DELDELREF":
                        return True
                    elif varsnv.triogenotype == "210" and varcnv.triogenotype == "DELREFDEL":
                        return True
                    else:
                        logging.info(
                            varid1 + " " + varid2 + " failed compound het screen, CNV and SNV pair "
                            "with triogenotypes not consistant with "
                            "compound het"
                        )
                else:
                    logging.info(
                        varid1 + " " + varid2 + " failed compound het screen, homozygous on "
                        "one alllle and other allele not CNV with "
                        "copy number of 1"
                    )

            elif (
                var1.get_mum_genotype() == "0"
                and var2.get_mum_genotype() != "0"
                and var1.get_dad_genotype() != "0"
                and var2.get_dad_genotype() == "0"
            ) or (
                var2.get_mum_genotype() == "0"
                and var1.get_mum_genotype() != "0"
                and var2.get_dad_genotype() != "0"
                and var1.get_dad_genotype() == "0"
            ):
                # one variant is inherited from each parent
                return True
            elif (
                (var1.dnm == True and var2.get_mum_genotype() != "0" and var2.get_dad_genotype() == "0")
                or (var1.dnm == True and var2.get_mum_genotype() == "0" and var2.get_dad_genotype() != "0")
                or (var2.dnm == True and var1.get_mum_genotype() != "0" and var1.get_dad_genotype() == "0")
                or (var2.dnm == True and var1.get_mum_genotype() == "0" and var1.get_dad_genotype() != "0")
            ):
                # one variant is DNM and the other is inherited
                return True
            else:
                logging.info(varid1 + " " + varid2 + " failed compound het screen")
                return False

        elif self.family.has_dad():
            pass
        # todo screening of compound hets with one parent
        elif self.family.has_mum():
            pass
        else:
            print("Error: should not get here as family must either have both, " "one or no parents")
            logging.error("Can't parse compound hets - family error")
            exit(1)
        logging.info(varid1 + " " + varid2 + " failed compound het screen")
        return False
