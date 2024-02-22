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

import json
from utils import params


def create_output(families, variants, inheritance_reports, outdir):
    """
    Create output file
    Identify variants in both compound het and single variants dicts
    Identify and flag possible MNVs
    Flag variants in cis in monoallelic genes
    """

    if len(families.keys()) > 1:
        outfile = outdir + "/" + "clinical_filter.tsv"
        inhreportfile = outdir + "/" + "clinical_filter_inheritance_report.txt"
    else:
        proband = families[list(families.keys())[0]].proband.person_id
        outfile = outdir + "/" + proband + "_clinical_filter.tsv"
        inhreportfile = outdir + "/" + proband + "_clinical_filter_inheritance_report.txt"

    header = [
        "family_id",
        "proband",
        "sex",
        "mum",
        "dad",
        "mum_aff",
        "dad_aff",
        "triogenotype",
        "decipher_genotype",
        "decipher_inheritance",
        "chrom",
        "pos",
        "stop",
        "ref",
        "alt",
        "DNM",
        "symbol",
        "result",
        "mode",
        "hgnc_id",
        "transcript",
        "canonical",
        "MANE_SELECT",
        "MANE_PLUS_CLINICAL",
        "consequence",
        "HGVSc",
        "HGVSp",
        "protein_position",
        "CLIN_SIG",
        "CALLSOURCE",
        "MEANLR2",
        "max_af",
        "ddd_af",
        "AC_het",
        "AC_hemi",
        "AC_tot",
        "GT",
        "GQ",
        "AD",
        "cnv_length",
        "cnv_copy_number",
        "mnv",
        "phased_15bp",
        "phased_any",
        "polyphen",
        "REVEL",
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
        "SpliceAI_pred",
        "LoF",
        "LoF_filter",
        "LoF_flags",
        "LoF_info",
        "CADD_PHRED",
    ]

    results = {}
    inhreports = {}

    for fam in variants.keys():
        phasedvars, phased_varids = identify_phased_variants(variants[fam])
        mnvs = identify_mnvs(phasedvars)
        variants_in_cis = identify_close_vars(phasedvars, 15)

        results[fam] = create_output_data(fam, families, variants, mnvs, variants_in_cis, phased_varids)

        inhreports[fam] = inheritance_reports[fam].__dict__

    print_output(results, header, outfile)
    print_inh_reports(inhreports, inhreportfile)


def print_output(results, header, outfile):
    with open(outfile, "w") as o:
        o.write(("\t").join(header))
        o.write("\n")

        for fam in results.keys():
            for var in results[fam].keys():
                line = ("\t").join(
                    [
                        results[fam][var]["family_id"],
                        results[fam][var]["proband"],
                        results[fam][var]["sex"],
                        results[fam][var]["mum"],
                        results[fam][var]["dad"],
                        results[fam][var]["mum_aff"],
                        results[fam][var]["dad_aff"],
                        results[fam][var]["triogenotype"],
                        results[fam][var]["decipher_genotype"],
                        results[fam][var]["decipher_inheritance"],
                        results[fam][var]["chrom"],
                        results[fam][var]["pos"],
                        results[fam][var]["stop"],
                        results[fam][var]["ref"],
                        results[fam][var]["alt"],
                        results[fam][var]["DNM"],
                        results[fam][var]["symbol"],
                        (",").join(results[fam][var]["result"]),
                        (",").join(results[fam][var]["mode"]),
                        results[fam][var]["hgnc_id"],
                        results[fam][var]["transcript"],
                        results[fam][var]["canonical"],
                        results[fam][var]["MANE_SELECT"],
                        results[fam][var]["MANE_PLUS_CLINICAL"],
                        results[fam][var]["consequence"],
                        results[fam][var]["HGVSc"],
                        results[fam][var]["HGVSp"],
                        results[fam][var]["protein_position"],
                        results[fam][var]["CLIN_SIG"],
                        results[fam][var]["CALLSOURCE"],
                        results[fam][var]["MEANLR2"],
                        results[fam][var]["max_af"],
                        results[fam][var]["ddd_af"],
                        results[fam][var]["AC_het"],
                        results[fam][var]["AC_hemi"],
                        results[fam][var]["AC_tot"],
                        results[fam][var]["GT"],
                        results[fam][var]["GQ"],
                        results[fam][var]["AD"],
                        results[fam][var]["cnv_length"],
                        results[fam][var]["cn"],
                        results[fam][var]["mnv"],
                        str(results[fam][var]["phased_15bp"]),
                        str(results[fam][var]["phased_any"]),
                        results[fam][var]["polyphen"],
                        results[fam][var]["REVEL"],
                        results[fam][var]["AlphaMissense_pred"],
                        results[fam][var]["AlphaMissense_rankscore"],
                        results[fam][var]["AlphaMissense_score"],
                        results[fam][var]["MPC_rankscore"],
                        results[fam][var]["MPC_score"],
                        results[fam][var]["PrimateAI_pred"],
                        results[fam][var]["PrimateAI_rankscore"],
                        results[fam][var]["PrimateAI_score"],
                        results[fam][var]["EVE_CLASS"],
                        results[fam][var]["EVE_SCORE"],
                        results[fam][var]["pLI_gene_value"],
                        results[fam][var]["SpliceAI_pred"],
                        results[fam][var]["LoF"],
                        results[fam][var]["LoF_filter"],
                        results[fam][var]["LoF_flags"],
                        results[fam][var]["LoF_info"],
                        results[fam][var]["CADD_PHRED"],
                    ]
                )
                o.write(line)
                o.write("\n")


def print_inh_reports(inhreports, inhreportfile):
    inhjson = json.dumps(inhreports)
    with open(inhreportfile, "w") as o:
        o.write(inhjson)
        o.write("\n")


def create_output_data(fam, families, variants, mnvs, variants_in_cis, phased_varids):
    # get family specific info
    results = {}

    if families[fam].has_mum():
        mum_aff = str(families[fam].mum.affected)
    else:
        mum_aff = "NA"
    if families[fam].has_dad():
        dad_aff = str(families[fam].dad.affected)
    else:
        dad_aff = "NA"

    # get variant specific info
    for varid in variants[fam]["single_variants"].keys():
        if varid in results.keys():
            results[varid]["result"].add("single_variant")
            results[varid]["mode"] = variants[fam]["single_variants"][varid]["mode"] | results[varid]["mode"]
        else:
            results[varid] = get_variant_info(
                variants[fam]["single_variants"][varid],
                varid,
                mnvs,
                variants_in_cis,
                phased_varids,
            )
            results[varid]["result"] = set(["single_variant"])
            fsplit = fam.split("_")
            results[varid]["family_id"] = fsplit[0]
            results[varid]["proband"] = fsplit[1]
            results[varid]["mum"] = families[fam].proband.mum_id
            results[varid]["dad"] = families[fam].proband.dad_id
            results[varid]["mum_aff"] = mum_aff
            results[varid]["dad_aff"] = dad_aff
            results[varid]["sex"] = families[fam].proband.sex

    for gn in variants[fam]["compound_hets"].keys():
        for varid in variants[fam]["compound_hets"][gn].keys():
            if varid in results.keys():
                results[varid]["result"].add("compound_het")
                results[varid]["mode"] = variants[fam]["compound_hets"][gn][varid]["mode"] | results[varid]["mode"]
            else:
                results[varid] = get_variant_info(
                    variants[fam]["compound_hets"][gn][varid],
                    varid,
                    mnvs,
                    variants_in_cis,
                    phased_varids,
                )
                results[varid]["result"] = set(["compound_het"])
                fsplit = fam.split("_")
                results[varid]["family_id"] = fsplit[0]
                results[varid]["proband"] = fsplit[1]
                results[varid]["mum"] = families[fam].proband.mum_id
                results[varid]["dad"] = families[fam].proband.dad_id
                results[varid]["mum_aff"] = mum_aff
                results[varid]["dad_aff"] = dad_aff
                results[varid]["sex"] = families[fam].proband.sex

    return results


def identify_phased_variants(variants):
    """
    Identify variants which are phased, these are then used to detect MNVs
    and variants in close proximity
    """
    phased = {}
    phased_varids = {}
    for v in variants["single_variants"].keys():
        if variants["single_variants"][v]["variant"].pid == ".":
            continue
        else:
            if not variants["single_variants"][v]["variant"].pid in phased.keys():
                phased[variants["single_variants"][v]["variant"].pid] = {}
            phased[variants["single_variants"][v]["variant"].pid][v] = {}
            phased[variants["single_variants"][v]["variant"].pid][v]["position"] = variants["single_variants"][v][
                "variant"
            ].pos
            phased[variants["single_variants"][v]["variant"].pid][v]["protein_position"] = variants["single_variants"][
                v
            ]["variant"].protein_position
            phased_varids[v] = 1

    for g in variants["compound_hets"].keys():
        for v in variants["compound_hets"][g].keys():
            if variants["compound_hets"][g][v]["variant"].pid == ".":
                continue
            else:
                if not variants["compound_hets"][g][v]["variant"].pid in phased.keys():
                    phased[variants["compound_hets"][g][v]["variant"].pid] = {}
                phased[variants["compound_hets"][g][v]["variant"].pid][v] = {}
                phased[variants["compound_hets"][g][v]["variant"].pid][v]["position"] = variants["compound_hets"][g][v][
                    "variant"
                ].pos
                phased[variants["compound_hets"][g][v]["variant"].pid][v]["protein_position"] = variants[
                    "compound_hets"
                ][g][v]["variant"].protein_position
                phased_varids[v] = 1

    return phased, phased_varids


def identify_close_vars(phasedvars, distance):
    """
    Identify variants in cis and in close proximity
    """
    close_vars = {}
    if len(phasedvars.keys()) >= 1:
        for pid in phasedvars.keys():
            if len(phasedvars[pid].keys()) > 1:
                for v1 in phasedvars[pid].keys():
                    pos1 = int(phasedvars[pid][v1]["position"])
                    for v2 in phasedvars[pid].keys():
                        pos2 = int(phasedvars[pid][v2]["position"])
                        if pos1 == pos2:
                            continue
                        else:
                            posdiff = abs(pos1 - pos2)
                            if posdiff <= distance:
                                close_vars[v1] = 1
                                close_vars[v2] = 1

    return close_vars


def identify_mnvs(phasedvars):
    """
    Identify MNVS which are in the same codon
    """
    mnvs = {}
    for vargroup in phasedvars.keys():
        if len(phasedvars[vargroup].keys()) >= 2:
            for v1 in phasedvars[vargroup].keys():
                pos1 = phasedvars[vargroup][v1]["position"]
                pp1 = phasedvars[vargroup][v1]["protein_position"]
                for v2 in phasedvars[vargroup].keys():
                    pos2 = phasedvars[vargroup][v2]["position"]
                    pp2 = phasedvars[vargroup][v2]["protein_position"]
                    if pos1 == pos2:
                        continue
                    else:
                        if pp1 == pp2:
                            mnvs[v1] = 1
                            mnvs[v2] = 1

    return mnvs


def get_variant_info(var, varid, mnvs, variants_in_cis, phased_varids):
    """
    Get variant specific information to go in the output lines
    """
    res = {}
    res["chrom"] = var["variant"].chrom
    res["pos"] = var["variant"].pos
    res["stop"] = stop_position(var)
    res["ref"] = var["variant"].ref
    res["alt"] = var["variant"].alt
    if var["variant"].is_snv():
        res["symbol"] = var["variant"].symbol
        res["hgnc_id"] = var["variant"].hgnc_id
    elif var["variant"].is_cnv():
        if len(var["variant"].reportable_symbol) == 0:
            res["symbol"] = "."
            res["hgnc_id"] = "."
        else:
            res["symbol"] = ("|").join(var["variant"].reportable_symbol)
            res["hgnc_id"] = ("|").join(var["variant"].reportable_hgnc_id)
    res["transcript"] = var["variant"].feature
    res["canonical"] = var["variant"].canonical
    res["MANE_SELECT"] = var["variant"].mane
    res["MANE_PLUS_CLINICAL"] = var["variant"].mane_clinical
    res["HGVSc"] = var["variant"].hgvsc
    res["HGVSp"] = var["variant"].hgvsp
    res["consequence"] = var["variant"].consequence
    res["protein_position"] = var["variant"].protein_position

    res["max_af"] = var["variant"].max_af
    res["ddd_af"] = var["variant"].ddd_af
    if var["variant"].is_snv():
        res["AC_het"] = var["variant"].AC_het
        res["AC_hemi"] = var["variant"].AC_hemi
        res["AC_tot"] = var["variant"].AC_tot
    else:
        res["AC_het"] = "."
        res["AC_hemi"] = "."
        res["AC_tot"] = "."
    res["GT"] = var["variant"].gt
    res["GQ"] = var["variant"].gq
    res["AD"] = var["variant"].ad
    res["cnv_length"] = var["variant"].cnv_length
    res["cn"] = var["variant"].cn
    res["triogenotype"] = var["variant"].triogenotype
    res["decipher_genotype"] = decipher_ready_genotype(var)
    res["decipher_inheritance"] = decipher_ready_inheritance(var)
    res["DNM"] = str(var["variant"].dnm)
    res["CLIN_SIG"] = var["variant"].CLIN_SIG
    res["CALLSOURCE"] = var["variant"].CALLSOURCE
    res["MEANLR2"] = var["variant"].MEANLR2

    if varid in mnvs.keys():
        res["mnv"] = "True"
    else:
        res["mnv"] = "False"
    if varid in variants_in_cis.keys():
        res["phased_15bp"] = "True"
    else:
        res["phased_15bp"] = "False"
    if varid in phased_varids.keys():
        res["phased_any"] = "True"
    else:
        res["phased_any"] = "False"

    res["mode"] = var["mode"]

    # CEPs
    res["polyphen"] = var["variant"].polyphen
    res["REVEL"] = var["variant"].revel
    res["AlphaMissense_pred"] = var["variant"].AlphaMissense_pred
    res["AlphaMissense_rankscore"] = var["variant"].AlphaMissense_rankscore
    res["AlphaMissense_score"] = var["variant"].AlphaMissense_score
    res["MPC_rankscore"] = var["variant"].MPC_rankscore
    res["MPC_score"] = var["variant"].MPC_score
    res["PrimateAI_pred"] = var["variant"].PrimateAI_pred
    res["PrimateAI_score"] = var["variant"].PrimateAI_score
    res["PrimateAI_rankscore"] = var["variant"].PrimateAI_rankscore
    res["EVE_CLASS"] = var["variant"].EVE_CLASS
    res["EVE_SCORE"] = var["variant"].EVE_SCORE
    res["pLI_gene_value"] = var["variant"].pLI_gene_value
    res["SpliceAI_pred"] = max_spliceAI(var)
    res["LoF"] = var["variant"].LoF
    res["LoF_filter"] = var["variant"].LoF_filter
    res["LoF_flags"] = var["variant"].LoF_flags
    res["LoF_info"] = var["variant"].LoF_info
    res["CADD_PHRED"] = var["variant"].CADD_PHRED

    return res


def decipher_ready_genotype(var):
    """
    Returns decipher ready genotype information for a given variant

    TODO :
    - CNV annotation
    - PAR regions annotation

    Args:
        var (dict): stores variant information

    Returns:
        str: decipher ready genotype information
    """

    triogenotype = var["variant"].triogenotype
    chrom = var["variant"].chrom
    sex = var["variant"].sex

    # CNVs are not handled
    decipher_genotype = params.DECIPHER_GENOTYPE_NA

    # Homozygous variants
    if triogenotype.startswith("2"):
        # For male with X chromosome variant, we might need to check on PAR region. See with Caroline.
        if (chrom == "X") and (sex == "XY"):
            decipher_genotype = params.DECIPHER_GENOTYPE_NA
        else:
            decipher_genotype = params.DECIPHER_GENOTYPE_HOMOZYGOUS

    # Heterozygous variants
    if triogenotype.startswith("1"):
        # For male with X chromosome variant, we might need to check on PAR region. See with Caroline.
        if (chrom == "X") and (sex == "XY"):
            decipher_genotype = params.DECIPHER_GENOTYPE_NA
        else:
            decipher_genotype = params.DECIPHER_GENOTYPE_HETEROZYGOUS

    return decipher_genotype


def decipher_ready_inheritance(var):
    """
    Returns decipher ready inheritance information for a given variant

    TODO :
    - CNV annotation

    Args:
        var (dict): stores variant information

    Returns:
        str: decipher ready inheritance information
    """

    triogenotype = var["variant"].triogenotype
    allelic_depths = var["variant"].ad

    # If the patient is a singleton, we do not have any clue about its inheritance
    if triogenotype.endswith("NANA"):
        return params.DECIPHER_INHERITANCE_UNKNOWN

    # CNVS are not handled at the moment
    if ("DUP" in triogenotype) or ("DEL" in triogenotype):
        return params.DECIPHER_INHERITANCE_NA

    # Denovo variants
    if triogenotype in ["100", "200"]:
        decipher_inheritance = params.DECIPHER_INHERITANCE_DENOVO

        if is_mosaic(allelic_depths):
            decipher_inheritance = params.DECIPHER_INHERITANCE_DENOVO_MOSAIC

    # Inherited variants
    else:

        # Biparental
        if triogenotype.endswith("11"):
            decipher_inheritance = params.DECIPHER_INHERITANCE_BIPARENTAL
        # Maternal
        elif triogenotype[1] != "0":
            decipher_inheritance = params.DECIPHER_INHERITANCE_MATERNAL

            if is_mosaic(allelic_depths):
                decipher_inheritance = params.DECIPHER_INHERITANCE_MATERNAL_MOSAIC
        # Paternal
        elif triogenotype[2] != "0":

            decipher_inheritance = params.DECIPHER_INHERITANCE_PATERNAL

            if is_mosaic(allelic_depths):
                decipher_inheritance = params.DECIPHER_INHERITANCE_PATERNAL_MOSAIC

    return decipher_inheritance


def is_mosaic(allelic_depths):
    """
    Given allelic depths for reference and alternalte allele,
    determine whether the variant could be mosaic

    Args:
        allelic_depths (tuple): reference and alternate alleles counts

    Returns:
        bool: mosaic or not
    """

    ref_allele_count = int(allelic_depths[0])
    alt_allele_count = int(allelic_depths[1])
    tot_allele_count = ref_allele_count + alt_allele_count

    return (alt_allele_count / tot_allele_count) < params.THRESHOLD_AD_MOSAICITY


def max_spliceAI(var):
    """
    SpliceAI output is quite lengthy (9 columns).
    This function retrieves the maximum spliceAI score across all evaluated splicing types
    and combine the results into a string (one column only)

    Args:
        var (dict): stores variant information

    Returns:
        str: splice_type:splice_score:splice_delta_position:splice_gene
    """

    splice_type = ["AG", "AL", "DG", "DL"]
    splice_delta_pos = [
        var["variant"].SpliceAI_pred_DP_AG,
        var["variant"].SpliceAI_pred_DP_AL,
        var["variant"].SpliceAI_pred_DP_DG,
        var["variant"].SpliceAI_pred_DP_DL,
    ]

    splice_scores = [
        var["variant"].SpliceAI_pred_DS_AG,
        var["variant"].SpliceAI_pred_DS_AL,
        var["variant"].SpliceAI_pred_DS_DG,
        var["variant"].SpliceAI_pred_DS_DL,
    ]
    splice_scores = [float(x) if x != "." else -1 for x in splice_scores]

    max_score, max_index = max((score, idx) for idx, score in enumerate(splice_scores))

    if max_score == -1:
        spliceAI_pred = "."
    else:
        spliceAI_pred = (
            f"{splice_type[max_index]}:{max_score}:{splice_delta_pos[max_index]}:{var['variant'].SpliceAI_pred_SYMBOL}"
        )

    return spliceAI_pred


def stop_position(var):
    """
    Determines stop position for CNVs and short deletions

    Args:
        var (dict): stores variant information

    Returns:
        str: stop position for CNVs and short deletions, "." otherwise
    """

    stop_pos = "."
    if var["variant"].is_cnv():
        stop_pos = int(var["variant"].pos) + int(var["variant"].cnv_length)
    else:
        length_ref = len(var["variant"].ref)
        if length_ref > 1:
            stop_pos = int(var["variant"].pos) + length_ref

    return str(stop_pos)
