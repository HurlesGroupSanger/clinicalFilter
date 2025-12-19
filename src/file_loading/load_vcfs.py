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

import subprocess
import logging
import os

from variants.snv import SNV
from variants.cnv import CNV


def load_variants(family, outdir, regions=None):
    """
    get variants in child and parents
    """
    proband_vcf = family.proband.get_vcf_path()
    if regions:
        regionfile = outdir + "/reg.tmp"
        with open(regionfile, "w") as rf:
            for r in regions:
                rf.write(r + "\n")

        sortedregs = outdir + "/reg.tmp_sorted"
        sortcmd = "sort -k1,1V -k2,2n -k3,3n " + regionfile + " > " + sortedregs
        os.system(sortcmd)
        child_vars = readvcf(proband_vcf, sortedregs, family.proband.get_sex())
        # remove regions files used for bcftools queries
        os.system("rm " + regionfile)
        os.system("rm " + sortedregs)
    else:
        child_vars = readvcf(proband_vcf, None, family.proband.get_sex())

    mum_vars = {}
    dad_vars = {}
    if not family.has_no_parents():
        # get a region string from the child variants as we don't need parental
        # variants which are not in the child

        childregionfile = outdir + "/childreg.tmp"
        with open(childregionfile, "w") as crf:
            for varid in child_vars.keys():
                idsplit = varid.split("_")
                reg = idsplit[0] + "\t" + idsplit[1]
                crf.write(reg + "\n")

        childsortedregs = outdir + "/childreg.tmp_sorted"
        childsortcmd = "sort -k1,1V -k2,2n " + childregionfile + " > " + childsortedregs

        os.system(childsortcmd)

        if family.has_mum():
            mum_vcf = family.mum.get_vcf_path()
            mum_vars = readvcf(mum_vcf, childsortedregs, "F")

        if family.has_dad():
            dad_vcf = family.dad.get_vcf_path()
            dad_vars = readvcf(dad_vcf, childsortedregs, "M")

        # remove child regions files used for bcftools queries
        os.system("rm " + childregionfile)
        os.system("rm " + childsortedregs)

    variants = {"child": child_vars, "mum": mum_vars, "dad": dad_vars}

    return variants


def readvcf(filename, regions, sex):
    """
    read vcf files and return a dict of variant objects
    """
    vars = {}

    # get list of info fields in the vcf
    infofields_wanted = [
        "Consequence",
        "Gene",
        "SYMBOL",
        "Feature",
        "CANONICAL",
        "MANE_SELECT",
        "MANE_PLUS_CLINICAL",
        "HGNC_ID",
        "MAX_AF",
        "MAX_AF_POPS",
        "DDD_AF",
        "DDD_father_AF",
        "REVEL",
        "PolyPhen",
        "Protein_position",
        "HGVSc",
        "HGVSp",
        "DNM",
        "DNG",
        "VAF",
        "END",
        "SVTYPE",
        "SVLEN",
        "CNVFILTER",
        "HGNC_ID_ALL",
        "SYMBOL_ALL",
        "AC_joint_XX",
        "AN_joint_XX",
        "nhomalt_joint_XX",
        "AC_joint_XY",
        "AN_joint_XY",
        "nhomalt_joint_XY",
        "AlphaMissense_pred",
        "AlphaMissense_rankscore",
        "AlphaMissense_score",
        "MPC_rankscore",
        "MPC_score",
        "PrimateAI_pred",
        "PrimateAI_rankscore",
        "PrimateAI_score",
        "EVE_CLASS",
        "EVE_SCORE",
        "pLI_gene_value",
        "SpliceAI_pred_DP_AG",
        "SpliceAI_pred_DP_AL",
        "SpliceAI_pred_DP_DG",
        "SpliceAI_pred_DP_DL",
        "SpliceAI_pred_DS_AG",
        "SpliceAI_pred_DS_AL",
        "SpliceAI_pred_DS_DG",
        "SpliceAI_pred_DS_DL",
        "SpliceAI_pred_SYMBOL",
        "LoF",
        "LoF_filter",
        "LoF_flags",
        "LoF_info",
        "CADD_PHRED",
        "CLIN_SIG",
        "ClinVar",
        "ClinVar_ALLELEID",
        "ClinVar_CLNSIG",
        "CALLSOURCE",
        "MEANLR2",
    ]
    formatfields = ["GT", "GQ", "PID", "AD", "CIFER_INHERITANCE", "CN"]

    # create infostring containing only the fields present
    info_query = []
    for inf in infofields_wanted:
        info_query.append("%INFO/" + inf)

    infostring = ("\t").join(info_query)
    formatstring = ("\t%").join(formatfields)

    if regions is None:
        bcfcmdroot = (
            "bcftools norm -m - "
            + filename
            + " | bcftools view -e 'INFO/MAX_AF>0.005 | FORMAT/GT[0]="
            + '"ref"'
            + "'  | bcftools query -u -f '%CHROM\t%POS\t%REF\t%ALT{0}\t"
        )
    else:
        bcfcmdroot = (
            "bcftools norm -m - -R "
            + regions
            + " "
            + filename
            + " | bcftools view -e 'INFO/MAX_AF>0.005 | FORMAT/GT[0]="
            + '"ref"'
            + "'  | bcftools query -u -f '%CHROM\t%POS\t%REF\t%ALT{0}\t"
        )

    bcfcmd = bcfcmdroot + infostring + "[\t%" + formatstring + "]\n'"
    output = runcommand(bcfcmd)

    if not output == "Error in command":
        outputlines = output.split("\n")
    else:
        logging.error("Variants not loaded from " + filename)

    for ol in outputlines:
        oldata = ol.split("\t")
        if len(oldata) < 2:
            continue
        alt = oldata[3]
        if alt == "*":  # get rid of any where alt allele is *
            continue
        # populate hash with variant data
        varid = ("_").join([oldata[0], oldata[1], oldata[2], alt])
        vdata = {}
        vdata["chrom"] = oldata[0]
        vdata["pos"] = oldata[1]
        vdata["ref"] = oldata[2]
        vdata["alt"] = alt
        vdata["consequence"] = oldata[4]
        vdata["ensg"] = oldata[5]
        vdata["symbol"] = oldata[6]
        vdata["feature"] = oldata[7]
        vdata["canonical"] = oldata[8]
        vdata["mane"] = oldata[9]
        vdata["mane_clinical"] = oldata[10]
        vdata["hgnc_id"] = oldata[11]
        vdata["max_af"] = oldata[12]
        vdata["max_af_pops"] = oldata[13]
        vdata["ddd_af"] = oldata[14]
        vdata["ddd_father_af"] = oldata[15]
        vdata["revel"] = oldata[16]
        vdata["polyphen"] = oldata[17]
        vdata["protein_position"] = oldata[18]
        vdata["hgvsc"] = oldata[19]
        vdata["hgvsp"] = oldata[20]
        vdata["sex"] = sex
        vdata["DNM"] = oldata[21]
        vdata["DNG"] = oldata[22]
        vdata["vaf"] = oldata[23]
        vdata["cnv_end"] = oldata[24]
        vdata["cnv_type"] = oldata[25]
        vdata["cnv_length"] = oldata[26]
        vdata["cnv_filter"] = oldata[27]
        vdata["hgnc_id_all"] = oldata[28]
        vdata["symbol_all"] = oldata[29]
        vdata["ac_XX"] = oldata[30]
        vdata["an_XX"] = oldata[31]
        vdata["nhomalt_XX"] = oldata[32]
        vdata["ac_XY"] = oldata[33]
        vdata["an_XY"] = oldata[34]
        vdata["nhomalt_XY"] = oldata[35]

        # Added for b38v3
        vdata["AlphaMissense_pred"] = oldata[36]
        vdata["AlphaMissense_rankscore"] = oldata[37]
        vdata["AlphaMissense_score"] = oldata[38]
        vdata["MPC_rankscore"] = oldata[39]
        vdata["MPC_score"] = oldata[40]
        vdata["PrimateAI_pred"] = oldata[41]
        vdata["PrimateAI_rankscore"] = oldata[42]
        vdata["PrimateAI_score"] = oldata[43]
        vdata["EVE_CLASS"] = oldata[44]
        vdata["EVE_SCORE"] = oldata[45]
        vdata["pLI_gene_value"] = oldata[46]
        vdata["SpliceAI_pred_DP_AG"] = oldata[47]
        vdata["SpliceAI_pred_DP_AL"] = oldata[48]
        vdata["SpliceAI_pred_DP_DG"] = oldata[49]
        vdata["SpliceAI_pred_DP_DL"] = oldata[50]
        vdata["SpliceAI_pred_DS_AG"] = oldata[51]
        vdata["SpliceAI_pred_DS_AL"] = oldata[52]
        vdata["SpliceAI_pred_DS_DG"] = oldata[53]
        vdata["SpliceAI_pred_DS_DL"] = oldata[54]
        vdata["SpliceAI_pred_SYMBOL"] = oldata[55]
        vdata["LoF"] = oldata[56]
        vdata["LoF_filter"] = oldata[57]
        vdata["LoF_flags"] = oldata[58]
        vdata["LoF_info"] = oldata[59]
        vdata["CADD_PHRED"] = oldata[60]
        vdata["CLIN_SIG"] = oldata[61]

        # Added for b38v4
        vdata["ClinVar"] = oldata[62]
        vdata["ClinVar_ALLELEID"] = oldata[63]
        vdata["ClinVar_CLNSIG"] = oldata[64]

        # Extra informations on CNVs
        vdata["CALLSOURCE"] = oldata[65]
        vdata["MEANLR2"] = oldata[66]

        # Format information
        vdata["gt"] = oldata[67]
        vdata["gq"] = oldata[68]
        vdata["pid"] = oldata[69]
        vdata["ad"] = oldata[70]
        vdata["cnv_inh"] = oldata[71]
        vdata["cn"] = oldata[72]

        if not vdata["DNM"] == "." or not vdata["DNG"] == ".":
            vdata["dnm"] = True
        else:
            vdata["dnm"] = False

        var = SNV
        if alt in ["<DEL>", "<DUP>"]:
            var = CNV
        if alt in ["<DEL>", "<DUP>"] and vdata["chrom"] == "Y":
            # exclude CNVs on Y
            logging.info(vdata["chrom"] + "_" + vdata["pos"] + "_" + vdata["ref"] + " CNV in Y: failed")
            continue
        vars[varid] = var(vdata)

    logging.info("Variants loaded from " + filename)

    return vars


def runcommand(cmd):
    try:
        byteoutput = subprocess.check_output(cmd, shell=True)
        return byteoutput.decode("UTF-8").rstrip()
    except subprocess.CalledProcessError:
        logging.debug(cmd)
        return "Error in command"
