#!/usr/bin/env python
import click
import pandas as pd
import textwrap
import os


@click.command()
@click.argument("list_vcfs")
@click.argument("list_probands")
@click.argument("families_ped")
@click.argument("cf_dir")
def create_ped(list_vcfs, list_probands, families_ped, cf_dir):
    """
    Create individual PED file used for clinical filtering

    Args:
        list_vcfs (str): list of VCFs file to consider in CF, output of pre-CF steps
        list_probands (str): list of probands to consider in CF, found in resources folder
        families_ped (str): relationships between individuals, found in resources folder
        cf_dir (str): clinical filtering output directory
    """

    # Load list of VCFs
    vcf_df = load_list_vcfs(list_vcfs)

    # Get probands to consider in CF
    probands_df = pd.read_csv(list_probands, header=None)
    probands_df.columns = ["stable_id"]
    probands_in_vcf_list = list(set(probands_df.stable_id) & set(vcf_df.stable_id))
    probands_df = probands_df.loc[probands_df.stable_id.isin(probands_in_vcf_list)]

    # Get relationships for the proband to consider
    family_df = load_families_ped(families_ped, probands_in_vcf_list)

    create_individual_ped(vcf_df, probands_in_vcf_list, family_df, cf_dir)


def load_list_vcfs(list_vcfs):
    """
    Load list of VCFs and perform some checking

    Args:
        list_vcfs (str): list of VCFs file to consider in CF
    """

    vcf_df = pd.read_csv(list_vcfs, header=None)
    vcf_df.columns = ["vcf_path"]
    vcf_df["stable_id"] = [x[-3] for x in vcf_df["vcf_path"].str.split("/")]
    vcf_df = vcf_df[["stable_id", "vcf_path"]]
    assert vcf_df.stable_id.str.startswith("DDDP").all()

    return vcf_df


def load_families_ped(families_ped, probands_in_vcf_list):
    """
    Load families relationships and perform some checking

    Args:
        families_ped (str): relationships between individuals
        probands_in_vcf_list (list): proband stable ids found in list of VCFs
    """

    family_df = pd.read_csv(families_ped, sep="\t")
    family_df.drop("path", inplace=True, axis=1)

    family_ids = list(family_df.loc[family_df.individual_id.isin(probands_in_vcf_list)].family_id)
    family_df = family_df.loc[family_df.family_id.isin(family_ids)].copy()

    family_df.replace({"sex": {"M": "XY", "F": "XX"}}, inplace=True)

    return family_df


def create_individual_ped(vcf_df, probands_in_vcf_list, family_df, cf_dir):
    """
    Create a ped file for all probands to consider in CF

    Args:
        vcf_df (pd.DataFrame): vcfs to consider in CF
        probands_in_vcf_list (list): proband stable ids found in list of VCFs
        family_df (pd.DataFrame): relationships between individuals
        cf_dir (str): clinical filtering folder
    """

    dict_ped = dict()
    for _, min_family_df in family_df.groupby("family_id"):

        # We check that we have all the VCF for the family
        assert min_family_df.individual_id.isin(vcf_df.stable_id).all()

        # Retrieve the proband in the family
        proband_id = list(set(min_family_df.individual_id) & set(probands_in_vcf_list))
        assert len(proband_id) == 1
        proband_id = proband_id[0]

        # Add VCF paths
        min_family_pluspath_df = min_family_df.merge(
            vcf_df[["stable_id", "vcf_path"]], left_on="individual_id", right_on="stable_id"
        )
        min_family_pluspath_df = min_family_pluspath_df[list(min_family_df.columns) + ["vcf_path"]]
        assert min_family_pluspath_df.shape[0] == min_family_df.shape[0]

        # Write ped file in the clinical filtering folder
        ped_file = write_individual_ped(min_family_pluspath_df, proband_id, cf_dir)
        dict_ped[proband_id] = ped_file

    with open("list_ped.tsv", "w") as f:
        for key, value in dict_ped.items():
            f.write(f"{key}\t{value}\n")


def write_individual_ped(min_family_df, proband_id, cf_dir):
    """
    Write individual ped in clinical filtering folder

    Args:
        min_family_df (pd.DataFrame): proband family relationship
        proband_id (str): proband identifier
        cf_dir (str): clinical filtering folder

    Returns:
        str: path of the individual ped file
    """

    proband_id_split = textwrap.wrap(proband_id, 2)
    ped_dir = f"{cf_dir}/{proband_id_split[0]}/{proband_id_split[1]}/{proband_id_split[2]}/{proband_id_split[3]}/{proband_id_split[4]}/{proband_id}/clinical_filter"
    os.makedirs(ped_dir, exist_ok=True)
    ped_file = f"{ped_dir}/{proband_id}.ped"

    min_family_df.to_csv(ped_file, sep="\t", header=False, index=False)

    # Returning realpath otherwise nextflow uses relative paths
    return os.path.realpath(ped_file)


if __name__ == "__main__":
    create_ped()
