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

from utils.utils import common_elements


class PreInheritanceFiltering(object):
    """
    Pre-inheritance filtering
    iterate through variants found in a proband
    remove any where GQ < 40 if autosome and any DDD_AF > 0.005
    remove any without functional consequence
    return a dict of variants per hgnc_id
    """

    def __init__(self, variants):
        self.variants = variants

    def preinheritance_filter(self):
        variants_per_gene = self.create_variants_per_gene()
        self.revel_filter(variants_per_gene)
        self.dnms_filter(variants_per_gene)
        self.X_maf_filter(variants_per_gene)
        return variants_per_gene

    def create_variants_per_gene(self):
        variants_per_gene = {}
        consequences = [
            "frameshift_variant",
            "missense_variant",
            "splice_donor_variant",
            "splice_acceptor_variant",
            "start_lost",
            "stop_gained",
            "protein_altering_variant",
            "transcript_ablation",
            "transcript_amplification",
            "inframe_insertion",
            "inframe_deletion",
            "stop_lost",
        ]

        for v in self.variants["child"].keys():
            # we only want SNVs in variants per gene
            if not self.variants["child"][v].is_snv():
                continue

            # fail if child GQ < 40
            if int(self.variants["child"][v].gq) < 40 and self.variants["child"][v].chrom not in ["X", "Y"]:
                logging.info(v + " failed low GQ: " + self.variants["child"][v].gq)
                continue

            # fail if DDD_AF > 0.005 (gnomAD AF variants above this threshold
            # are not loaded)
            if float(self.variants["child"][v].ddd_af) > 0.005:
                logging.info(v + " failed high DDD AF: " + self.variants["child"][v].ddd_af)
                continue

            # Check if the variant has a high spliceAI score
            is_high_spliceAI_DNM = self.is_high_spliceAI_DNM(v)

            # If the variant is not a DNM with high spliceAI score, it has to have a functional consequence to be kept
            if not is_high_spliceAI_DNM:
                cqs = self.variants["child"][v].consequence.split("&")
                coding_cqs = common_elements(cqs, consequences)
                if len(coding_cqs) == 0:
                    logging.info(v + " failed, no functional consequences: " + self.variants["child"][v].consequence)
                    continue

            hgncid = cqs = self.variants["child"][v].hgnc_id
            if not hgncid in variants_per_gene.keys():
                variants_per_gene[hgncid] = {}

            variants_per_gene[hgncid][v] = {}
            variants_per_gene[hgncid][v]["child"] = self.variants["child"][v]
            if v in self.variants["mum"].keys():
                variants_per_gene[hgncid][v]["mum"] = self.variants["mum"][v]
            if v in self.variants["dad"].keys():
                variants_per_gene[hgncid][v]["dad"] = self.variants["dad"][v]

        return variants_per_gene

    def is_high_spliceAI_DNM(self, v):
        """
        Keep UTR and synonymous DNM variants if they have a high spliceAI score

        Args:
            v (str): variant identifier

        Returns:
            bool: whether or not the variant is a DNM with a high spliceAI score
        """

        if not self.variants["child"][v].dnm:
            return False

        SPLICE_AI_THRESHOLD = 0.8

        splice_ai_consequences = [
            "synonymous_variant",
            "3_prime_UTR_variant",
            "5_prime_UTR_variant",
        ]
        cqs = self.variants["child"][v].consequence.split("&")
        coding_cqs = common_elements(cqs, splice_ai_consequences)

        if coding_cqs:

            ds_ag = (
                float(self.variants["child"][v].SpliceAI_pred_DS_AG)
                if self.variants["child"][v].SpliceAI_pred_DS_AG != "."
                else 0
            )

            ds_al = (
                float(self.variants["child"][v].SpliceAI_pred_DS_AL)
                if self.variants["child"][v].SpliceAI_pred_DS_AL != "."
                else 0
            )

            ds_dg = (
                float(self.variants["child"][v].SpliceAI_pred_DS_DG)
                if self.variants["child"][v].SpliceAI_pred_DS_DG != "."
                else 0
            )

            ds_dl = (
                float(self.variants["child"][v].SpliceAI_pred_DS_DL)
                if self.variants["child"][v].SpliceAI_pred_DS_DL != "."
                else 0
            )

            if (
                (ds_ag >= SPLICE_AI_THRESHOLD)
                | (ds_al >= SPLICE_AI_THRESHOLD)
                | (ds_dg >= SPLICE_AI_THRESHOLD)
                | (ds_dl >= SPLICE_AI_THRESHOLD)
            ):
                return True

        return False

    def revel_filter(self, variants_per_gene):
        """
        Remove missense variants with REVEL < 0.5 unless DNM
        """
        for gn in list(variants_per_gene.keys()):
            for varid in list(variants_per_gene[gn].keys()):
                childvar = variants_per_gene[gn][varid]["child"]
                if childvar.dnm == True:
                    continue
                elif childvar.consequence.find("missense_variant") == -1:
                    continue
                elif childvar.revel == ".":
                    continue
                else:
                    revel = float(childvar.revel)
                    if revel < 0.4:
                        logging.info(varid + " failed REVEL filter: " + str(revel))
                        del variants_per_gene[gn][varid]
                        if len(variants_per_gene[gn].keys()) < 1:
                            del variants_per_gene[gn]

    def dnms_filter(self, variants_per_gene):
        """
        Remove DNMs that don't pass filters
        """
        for gn in list(variants_per_gene.keys()):
            for varid in list(variants_per_gene[gn].keys()):
                childvar = variants_per_gene[gn][varid]["child"]
                if (childvar.triogenotype == "100" or childvar.triogenotype == "200") and childvar.dnm == False:
                    logging.info(varid + " triogenotype = " + childvar.triogenotype + " and failed DNM filter")
                    del variants_per_gene[gn][varid]
                    if len(variants_per_gene[gn].keys()) < 1:
                        del variants_per_gene[gn]

    def X_maf_filter(self, variants_per_gene):
        """
        Variants in X have more stringent allele frequencies - fail if
        gnomad > 0.000001 or DDD unaffected father > 0
        This is in a separate function for clarity and ease of modification
        """
        for gn in list(variants_per_gene.keys()):
            for varid in list(variants_per_gene[gn].keys()):
                childvar = variants_per_gene[gn][varid]["child"]
                if childvar.chrom == "X":
                    max_af = childvar.max_af
                    if max_af == ".":
                        max_af = "0"
                    ddd_father_af = childvar.ddd_father_af
                    if ddd_father_af == ".":
                        ddd_father_af = "0"
                    if float(max_af) > 0.000001:
                        logging.info(varid + " failed X chromosome allele " "frequency: gnomad AF = " + str(max_af))
                        del variants_per_gene[gn][varid]
                        if len(variants_per_gene[gn].keys()) < 1:
                            del variants_per_gene[gn]
                    elif float(ddd_father_af) > 0:
                        logging.info(
                            varid + " failed X chromosome allele "
                            "frequency: DDD unaffected father "
                            "AF = " + str(ddd_father_af)
                        )
                        del variants_per_gene[gn][varid]
                        if len(variants_per_gene[gn].keys()) < 1:
                            del variants_per_gene[gn]
