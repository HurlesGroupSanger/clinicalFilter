# tests/helpers/ped_factory.py
import pandas as pd
from typing import Iterable, Tuple


def make_ped_text(family_id: str, members: Iterable[Tuple[str, str, str, int, int]]) -> str:
    """Return a PED text (no header). members: (id, father, mother, sex, phenotype)"""
    rows = [f"{family_id}\t{ind}\t{dad}\t{mom}\t{sex}\t{phen}" for ind, dad, mom, sex, phen in members]
    return "\n".join(rows) + "\n"


def add_absolute_paths(ped_file, base_path, out_ped):
    """
    Add absolute path to the VCF in the test PED file and create a new PED file.

    Args:
        ped_file (str): test PED file
        base_path (str): installation directory
        out_ped (str) : test PED file with absolute path
    """

    ped_df = pd.read_csv(ped_file, sep="\t", header=None)
    ped_df.iloc[:, 6] = ped_df.iloc[:, 6].apply(lambda x: f"{base_path}/{x}")

    ped_df.to_csv(out_ped, sep="\t", header=False, index=False)
