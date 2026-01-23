#!/usr/bin/env python

import click
import pandas as pd
import json
import logging
import logging.config


@click.command()
@click.argument("config_file")
def annotate_cf(config_file):
    """
    Annotate CF results

    Args:
        config_file (str): path to configuration file
    """

    logger = logging.getLogger("logger")

    conf = get_conf(config_file)

    # Load all files and construct variant identifiers
    logger.info("Loading resources...")
    id_mapping_df = load_id_mapping(conf["id_mapping"])
    previous_gene_list = load_previous_genelist(conf["previous_gene_list"])
    decipher_variants_info_df = load_decipher_variants_info(conf["decipher_variants_info"], id_mapping_df)
    if "b37_cf_results" in conf:
        b37_cf_results = load_b37_cf_results(conf["b37_cf_results"])
    else:
        b37_cf_results = pd.DataFrame()
    latest_cf_results = load_latest_cf_results(conf["latest_cf_results"])
    tiering_info = load_tiering_info(conf["tiering_info"])

    if "b38_cf_previous_results" in conf:
        b38_cf_previous_results = load_latest_cf_results(conf["b38_cf_previous_results"])
    else:
        b38_cf_previous_results = pd.DataFrame()

    # Annotate CF results
    annotated_df = annotate(
        latest_cf_results,
        b37_cf_results,
        decipher_variants_info_df,
        previous_gene_list,
        id_mapping_df,
        b38_cf_previous_results,
        tiering_info,
    )

    # Order columns, sort rows
    annotated_df = format_results(annotated_df)

    # Export results
    annotated_df.to_csv(
        f"{conf['outdir']}/{conf['latest_cf_results'].split('/')[-1].replace('.tsv', '_annotated.tsv')}",
        sep="\t",
        index=False,
    )

    # Export only new variants
    new_variants_df = keep_new_variants_only(annotated_df)
    new_variants_df.to_csv(
        f"{conf['outdir']}/{conf['latest_cf_results'].split('/')[-1].replace('.tsv', '_annotated_new_variants.tsv')}",
        sep="\t",
        index=False,
    )


def annotate(
    cf_results,
    b37_cf_results,
    decipher_variants_info,
    previous_gene_list,
    id_mapping_df,
    b38_cf_previous_results,
    tiering_info,
):
    """
    Annotate CF results with various information

    Args :
        cf_results (pd.DataFrame): current CF results
        b37_cf_results (pd.DataFrame): last B37 CF results
        decipher_variants_info (pd.DataFrame): variants already reported in DECIPHER
        previous_gene_list (list): DDG2P genes used in the previous run
        id_mapping_df (pd.DataFrame) : DDD identifiers mapping
        b38_cf_previous_results (pd.DataFrame): previous B38 CF results
        tiering_info (pd.DataFrame): diagnostic tiering information for each patient

    """

    logger = logging.getLogger("logger")
    logger.info("Starting CF results annotation...")

    cf_results["ref_reads"] = "."
    cf_results["alt_reads"] = "."
    cf_results["indel_length"] = "."
    cf_results["in_previous_build_38"] = "."
    if not b37_cf_results.empty:
        cf_results["in_build_37"] = "n"
    cf_results["in_decipher"] = "n"
    cf_results["vars_per_gene"] = 0
    cf_results["gene_in_prev_run"] = "n"
    cf_results["cnvs_per_proband"] = 0
    cf_results["cnv_chromosome_count"] = 0
    cf_results["vars_per_proband"] = 0

    # Add DECIPHER ids
    cf_results = cf_results.merge(
        id_mapping_df[["decipher_id", "person_stable_id"]],
        left_on="proband",
        right_on="person_stable_id",
    )
    cf_results.drop("person_stable_id", axis=1, inplace=True)

    # Compute indel length
    cf_results["indel_length"] = cf_results.apply(get_indel_length, axis=1)

    # Get allelic depths
    cf_results[["ref_reads", "alt_reads"]] = cf_results.apply(get_allelic_depth, axis=1)

    # Get number of variants in the list associated to the current proband
    nb_variants_per_proband = dict(cf_results["proband"].value_counts())
    cf_results["vars_per_proband"] = cf_results.apply(
        lambda row, nb_variants_per_proband: nb_variants_per_proband[row.proband],
        axis=1,
        args=(nb_variants_per_proband,),
    )

    # Check if another variant is present in the same region (+-100bp)
    cf_results = other_variant_in_same_region(cf_results)

    # Get number of variants in the list associated to the current gene
    nb_variants_per_gene = dict(cf_results["symbol"].value_counts())
    cf_results["vars_per_gene"] = cf_results.apply(
        lambda row, nb_variants_per_gene: nb_variants_per_gene[row.symbol],
        axis=1,
        args=(nb_variants_per_gene,),
    )

    probands = cf_results["proband"].unique()

    # Get number chromosome harborin a CNV for current proband
    nb_chrom_with_a_cnv_in_proband = dict()
    for proband in probands:
        nb_chrom_with_a_cnv_in_proband[proband] = cf_results.loc[
            (cf_results["proband"] == proband) & (cf_results["cnv"] == "y"), "chrom"
        ].nunique()
    cf_results["cnv_chromosome_count"] = cf_results.apply(
        lambda row, nb_chrom_with_a_cnv_in_proband: nb_chrom_with_a_cnv_in_proband[row.proband],
        axis=1,
        args=(nb_chrom_with_a_cnv_in_proband,),
    )

    # Get number chromosome harborin a CNV for current proband
    nb_cnv_in_proband = dict()
    for proband in probands:
        nb_cnv_in_proband[proband] = ((cf_results["proband"] == proband) & (cf_results["cnv"] == "y")).sum()
    cf_results["cnvs_per_proband"] = cf_results.apply(
        lambda row, nb_cnv_in_proband: nb_cnv_in_proband[row.proband],
        axis=1,
        args=(nb_cnv_in_proband,),
    )

    # Get number of CNVs in the proband
    nb_cnvs_per_proband = dict(zip(probands, [0] * len(probands)))
    for idx, row in cf_results.iterrows():
        if row.cnv == "y":
            nb_cnvs_per_proband[row.proband] += 1

    cf_results["cnvs_per_proband"] = cf_results["proband"].map(nb_cnvs_per_proband)

    # Check if variants was already in B37 results
    logger.info("Checking whether variants were already in previous CF results and/or in DECIPHER...")

    if not b37_cf_results.empty:
        cf_results["in_build_37"] = cf_results.apply(
            lambda row, b37_cf_results: ("y" if row.varid in list(b37_cf_results.varid) else "n"),
            axis=1,
            args=(b37_cf_results,),
        )

    # Check if variants was already in previous b38 results
    if not b38_cf_previous_results.empty:
        cf_results["in_previous_build_38"] = cf_results.apply(
            lambda row, b38_cf_previous_results: ("y" if row.varid in list(b38_cf_previous_results.varid) else "n"),
            axis=1,
            args=(b38_cf_previous_results,),
        )

    # Check if variants has been reported in DECIPHER
    cf_results["in_decipher"] = cf_results.apply(
        lambda row, decipher_variants_info: ("y" if row.varid in list(decipher_variants_info.varid) else "n"),
        axis=1,
        args=(decipher_variants_info,),
    )

    # Check if the variant's gene was already in the previous DDG2P list
    cf_results["gene_in_prev_run"] = cf_results.apply(
        lambda row, previous_gene_list: ("y" if row.symbol in list(previous_gene_list) else "n"),
        axis=1,
        args=(previous_gene_list,),
    )

    # Fuzzy matching for CNVs
    logger.info("CNV fuzzy matching...")
    if not b37_cf_results.empty:
        cf_results = cnv_fuzzy_matching(cf_results, b37_cf_results, "in_build_37")
    cf_results = cnv_fuzzy_matching(cf_results, decipher_variants_info, "in_decipher")
    cf_results = cnv_fuzzy_matching(cf_results, b38_cf_previous_results, "in_previous_build_38")

    # Fuzzy matching for indels/MNVs
    logger.info("Indels/MNVs fuzzy matching...")
    cf_results = indels_mnvs_fuzzy_matching(cf_results, decipher_variants_info)

    # Add diagnostic tiering information
    cf_results = cf_results.merge(tiering_info, on="decipher_id", how="left")

    cf_results.drop(["family_id", "cnv"], inplace=True, axis=1)

    return cf_results


def get_indel_length(row):
    if row.alt not in ["<DEL>", "<DUP>"]:
        lendiff = abs(len(row.ref) - len(row.alt))
        if not lendiff == 0:
            return str(lendiff)
    return "."


def get_allelic_depth(row):
    allelic_depths = row.AD.split(",")
    # Compound het variants have been split during CF so max number of alleles is 2
    if len(allelic_depths) > 1:
        allelic_depths = (allelic_depths[0], allelic_depths[-1])
    else:
        allelic_depths = (".", ".")

    return pd.Series({"ref_reads": allelic_depths[0], "alt_reads": allelic_depths[1]})


def get_nb_variants_per_proband(row, df):
    """_summary_

    Args:
        row (_type_): _description_
        df (_type_): _description_
    """

    return (df.proband == row.proband).sum()


def keep_new_variants_only(df):
    """
    Keep only new variants (not found in B37 or DECIPHER and not in previous B38 build).
    Exception for DNMs which are kept even if already selected in previous CF run but not reported in DECIPHER.

    Args:
        df (pd.DataFrame): CF results
    """

    def filter_cnvs(row):
        filter = True
        if row.alt in ["<DEL>", "<DUP>"]:
            if row.cnvs_per_proband > 4:
                filter = False
        return filter

    # Keep only variants that are not already in DECIPHER
    df = df.loc[df.in_decipher == "n"]

    # Keep only variants that were not already selected in last B37 run (or are de novo not reported in DECIPHER)
    if "in_build_37" in df.columns:
        filt = (df.in_build_37 == "n") | (df.decipher_inheritance.str.startswith("de_novo"))
        df = df.loc[filt]

    # Keep only variants that were not already selected in previous B38 build (or are de novo not reported in DECIPHER)
    if "in_previous_build_38" in df.columns:
        filt = (df.in_previous_build_38 == "n") | (df.decipher_inheritance.str.startswith("de_novo"))
        df = df.loc[filt]

    # Filter CNVS from probands having more than 4 CNVs
    df = df[df.apply(filter_cnvs, axis=1)]

    return df


def get_conf(config_file):
    """
    Load parameters from json file

    Args:
        config_file (str): path to configuration file

    Returns:
        dict: parameters
    """

    with open(config_file, "r") as f:
        conf = json.load(f)

    return conf


def load_id_mapping(filename):
    """
    Load the file that maps all DDD identifiers together

    Args:
        filename (str): path to mapping file

    Returns:
        pd.DataFrame: identifiers mapping
    """
    df = pd.read_csv(filename, sep="\t")
    return df


def load_previous_genelist(filename):
    """
    Load DDG2P genes used in the previous run

    Args:
        filename (str): path to the previous gene list
    """

    df = pd.read_csv(filename, sep="\t")

    # If the user provided a gene list
    if df.shape[1] < 2:
        df = pd.read_csv(filename, sep="\t", header=None)
        genelist = list(df.iloc[:, 0])
    # Or if he provided the DDG2P file
    else:
        genelist = list(df.gene.unique())

    return genelist


def load_decipher_variants_info(filename, id_mapping_df):
    """
    Load variants already reported in DECIPHER

    Args:
        filename (str): path to the DECIPHER variants file
        id_mapping_df (pd.DataFrame) : DDD identifiers mapping
    """

    df = pd.read_csv(filename, sep="\t", dtype={"patient_id": str}, low_memory=False)
    df = df.merge(id_mapping_df, left_on="patient_id", right_on="decipher_id")

    df = build_decipher_variant_id(df)
    return df


def build_decipher_variant_id(df):
    """
    Build DECIPHER variant ids by combining ther genomic coordinates, ref and alt alleles.

    Args:
        df (pd.DataFrame): DECIPHER variants

    Returns:
        pd.DataFrame : DECIPHER variants with built-in variant identifier
    """

    list_var_id = list()
    for idx, row in df.iterrows():
        variant_class = row.variant_class
        if variant_class == "sequence_variant":
            varid = ("_").join(
                [
                    row.person_stable_id,
                    str(row.chr),
                    str(row.start),
                    row.ref_allele,
                    row.alt_allele,
                ]
            )
        elif variant_class == "deletion":
            varid = ("_").join([row.person_stable_id, str(row.chr), str(row.start), "DEL"])
        elif variant_class == "duplication":
            varid = ("_").join([row.person_stable_id, str(row.chr), str(row.start), "DUP"])
        elif variant_class == "amplification":
            varid = ("_").join([row.person_stable_id, str(row.chr), str(row.start), "DUP"])
        else:
            varid = ""

        list_var_id.append(varid)

    df["varid"] = list_var_id

    return df


def load_b37_cf_results(filename):
    """
    Load last B37 CF results

    Args:
        filename (str): path to last B37 CF results
    """

    df = pd.read_csv(filename, sep="\t")
    df = build_b37_variant_id(df)
    return df


def build_b37_variant_id(df):
    """
    Build variant ids by combining ther genomic coordinates, ref and alt alleles.

    Args:
        df (pd.DataFrame): B37 variants

    Returns:
        pd.DataFrame: B37 variants with built-in variant identifiers
    """

    list_var_id = list()
    for idx, row in df.iterrows():
        ref = row["ref/alt_alleles"].split("/")[0]
        alt = row["ref/alt_alleles"].split("/")[1]
        if alt == "<DEL>":
            varid = ("_").join([row["proband"], str(row.chrom), str(row.position), "DEL"])
        elif alt == "<DUP>":
            varid = ("_").join([row["proband"], str(row.chrom), str(row.position), "DUP"])
        else:
            npos, nref, nalt = normalise_variant(row.position, ref, alt)
            varid = ("_").join([row["proband"], str(row.chrom), str(npos), nref, nalt])

        list_var_id.append(varid)

    df["varid"] = list_var_id

    return df


def normalise_variant(pos, ref, alt):
    """
    Normalise (left justify) a variant to aid variant matching


    Args:
        pos (int): variant genomic position
        ref (str): variant reference allele
        alt (str): variant alternate allele

    Returns:
        tuple : normalized pos, ref and alt
    """
    pos = int(pos)
    # If it's a simple SNV, don't remap anything
    if len(ref) == 1 and len(alt) == 1:
        return pos, ref, alt
    else:
        # strip off identical suffixes
        while alt[-1] == ref[-1] and min(len(alt), len(ref)) > 1:
            alt = alt[:-1]
            ref = ref[:-1]
        # strip off identical prefixes and increment position
        while alt[0] == ref[0] and min(len(alt), len(ref)) > 1:
            alt = alt[1:]
            ref = ref[1:]
            pos += 1
        return pos, ref, alt


def load_latest_cf_results(filename):
    """
    Load last B38 CF results

    Args:
        filename (str): path to last B37 CF results
    """

    df = pd.read_csv(filename, sep="\t", low_memory=False)
    df = build_b38_variant_id(df)
    return df


def build_b38_variant_id(df):
    """
    Build variant ids by combining ther genomic coordinates, ref and alt alleles.

    Args:
        df (pd.DataFrame): B38 variants
    """

    list_var_id = list()
    list_cnv = list()
    for idx, row in df.iterrows():
        if row.alt == "<DEL>":
            varid = ("_").join([row["proband"], str(row.chrom), str(row.pos), "DEL"])
            list_cnv.append("y")
        elif row.alt == "<DUP>":
            varid = ("_").join([row["proband"], str(row.chrom), str(row.pos), "DUP"])
            list_cnv.append("y")
        else:
            npos, nref, nalt = normalise_variant(row.pos, row.ref, row.alt)
            varid = ("_").join([row["proband"], str(row.chrom), str(npos), nref, nalt])
            list_cnv.append("n")

        list_var_id.append(varid)

    df["varid"] = list_var_id
    df["cnv"] = list_cnv

    return df


def load_tiering_info(filename):
    """
    Load diagnostic tiering information for each patient

    Args:
        filename (str): path to tiering information file
    """

    df = pd.read_csv(filename, sep="\t", low_memory=False, dtype={"decipher_id": str})
    return df


def other_variant_in_same_region(cf_df):
    """
    Flag variants occurring in the same region (+- 100bp) as another variant in the same proband

    Args:
        cf_df (pd.DataFrame): CF results
    """

    OFFSET = 100

    cf_df["other_variant_in_region"] = "n"
    for _, patient_variants_df in cf_df.groupby("decipher_id"):

        for idx, row in patient_variants_df.iterrows():
            chrom = row["chrom"]
            pos = row["pos"]

            overlapping_variants = patient_variants_df.loc[
                (patient_variants_df["chrom"] == chrom)
                & (patient_variants_df["pos"] >= (pos - OFFSET))
                & (patient_variants_df["pos"] <= (pos + OFFSET))
                & (patient_variants_df.index != idx)
            ]

            if not overlapping_variants.empty:
                cf_df.at[idx, "other_variant_in_region"] = "y"

    return cf_df


def cnv_fuzzy_matching(cf_df, other_df, column_name):
    """
    Check whether the CNV was previously reported with slightly different positions (using an offset of 1000bp)
    NOTE : Right now we are using CNV lifted over from b37 so no difference should be observed between
    last b38 and previous b38/b37 results, which won't be the case if we call CNVs on b38 directly

    Args:
        cf_results (pd.DataFrame): current run CF results
        other_df (pd.DataFrame): previous run CF results or variants reported in DECIPHER
        column_name (str) : column to use to flag CNVs already found/reported according to the fuzzy matching
    """

    OFFSET_CNV = 1000

    cnv_cf_df = cf_df.loc[cf_df.alt.isin(["<DEL>", "<DUP>"])]
    if column_name == "in_decipher":
        cnv_other_df = other_df.loc[other_df.variant_class.isin(["deletion", "duplication"])]
        proband_column = "person_stable_id"
        chrom_column = "chr"
    elif column_name == "in_build_37":
        cnv_other_df = other_df.loc[other_df.varid.str.contains("DEL") | other_df.varid.str.contains("DUP")]
        proband_column = "proband"
        chrom_column = "chrom"
    else:
        cnv_other_df = other_df.loc[other_df.alt.isin(["<DEL>", "<DUP>"])]
        proband_column = "proband"
        chrom_column = "chrom"

    idx_cnv_match = list()
    for idx, row in cnv_cf_df.iterrows():
        proband_id = row["proband"]
        for idx2, row2 in cnv_other_df.loc[cnv_other_df[proband_column] == proband_id].iterrows():
            if str(row2[chrom_column]) != str(row.chrom):
                continue
            if column_name == "in_decipher":
                if (abs(row["pos"] - row2["start"]) < OFFSET_CNV) | (
                    abs(row["pos"] + int(row["cnv_length"]) - row2["end"]) < OFFSET_CNV
                ):
                    idx_cnv_match.append(idx)
            elif column_name == "in_build_37":
                if (abs(row["pos"] - row2["position"]) < OFFSET_CNV) | (
                    abs(row["pos"] + int(row["cnv_length"]) - (row2["position"] + row2["cnv_length"])) < OFFSET_CNV
                ):
                    idx_cnv_match.append(idx)
            else:
                if (abs(row["pos"] - row2["pos"]) < OFFSET_CNV) | (
                    abs(row["pos"] + int(row["cnv_length"]) - (row2["pos"] + int(row2["cnv_length"]))) < OFFSET_CNV
                ):
                    idx_cnv_match.append(idx)

    cf_df.loc[idx_cnv_match, column_name] = "y"

    return cf_df


def indels_mnvs_fuzzy_matching(cf_df, decipher_variants_info):
    """
    Check whether indels/MNVs were previously reported with slightly different positions

    Args:
        cf_results (pd.DataFrame): current run CF results
    """

    logger = logging.getLogger("logger")

    OFFSET_INDEL_MNV = 20

    # Exclude CNVs
    sequence_variants_df = cf_df.loc[~cf_df.alt.isin(["<DEL>", "<DUP>"])]

    patients_with_variant_not_in_decipher = sequence_variants_df.loc[
        sequence_variants_df.in_decipher == "n", "decipher_id"
    ].unique()

    # For each proband, check if any variant match with DECIPHER reported variants
    for decipher_id in patients_with_variant_not_in_decipher:

        patient_variants_cf_df = sequence_variants_df.loc[sequence_variants_df.decipher_id == decipher_id]
        patient_decipher_variants_df = decipher_variants_info.loc[decipher_variants_info.decipher_id == decipher_id]

        # If no overlapping chromosomes between variants already reported in DECIPHER and the ones not yet reported in DECIPHER, skip
        if (
            set(patient_decipher_variants_df.chr)
            & set(patient_variants_cf_df.loc[patient_variants_cf_df.in_decipher == "n"].chrom)
        ) == set():
            continue

        # For each variant not yet reported in DECIPHER, check if any variant reported in DECIPHER is close enough
        for idx, patient_variant in patient_variants_cf_df.loc[patient_variants_cf_df.in_decipher == "n"].iterrows():

            patient_decipher_variants_df_on_chrom = patient_decipher_variants_df.loc[
                patient_decipher_variants_df.chr == patient_variant.chrom
            ]

            for _, decipher_variant in patient_decipher_variants_df_on_chrom.iterrows():
                if abs(patient_variant["pos"] - decipher_variant["start"]) < OFFSET_INDEL_MNV:
                    cf_df.at[idx, "in_decipher"] = "y"
                    logger.info("Found fuzzy match for {}".format(patient_variant.varid))
                    break

    return cf_df


def format_results(df):
    """
    Format CF results

    Args:
        df (pd.DataFrame): CF results
    """

    # Sort results
    df.sort_values(["decipher_id", "chrom", "pos"], inplace=True)

    # Remove allelic depth column as we created a ref_reads and alt_reads columns
    df.drop("AD", axis=1, inplace=True)

    # Order columns
    after_column = "LoF_info"
    insert_after_column = "phased_any"

    after_column_idx = df.columns.get_loc(after_column)
    insert_after_column_idx = df.columns.get_loc(insert_after_column)

    main_columns = list(df.columns[: insert_after_column_idx + 1])
    ceps_columns = list(df.columns[insert_after_column_idx + 1 : after_column_idx + 1])
    postcf_columns = list(df.columns[after_column_idx + 1 :])

    df = df[main_columns + postcf_columns + ceps_columns]

    # Put decipher identifier in front
    df = df[["decipher_id"] + [x for x in df.columns if x not in ["decipher_id", "varid"]] + ["varid"]]

    return df


def init_log(show_time=True):
    """
    Initialise logger

    Args:
        show_time (bool, optional): print time in log entries. Defaults to True.
    """

    if show_time:
        log_format = "[%(levelname)s:%(asctime)s] %(message)s"
    else:
        log_format = "%(levelname)s | %(message)s"

    MY_LOGGING_CONFIG = {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "default_formatter": {"format": log_format},
        },
        "handlers": {
            "stream_handler": {
                "class": "logging.StreamHandler",
                "formatter": "default_formatter",
            },
        },
        "loggers": {
            "logger": {
                "handlers": ["stream_handler"],
                "level": "INFO",
                "propagate": True,
            }
        },
    }

    logging.config.dictConfig(MY_LOGGING_CONFIG)


if __name__ == "__main__":
    init_log()
    annotate_cf()
