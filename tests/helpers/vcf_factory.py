# tests/helpers/vcf_factory.py
from textwrap import dedent
from typing import Iterable, Tuple


VCF_HEADER = dedent(
    """\
##fileformat=VCFv4.2
##source=test
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""
)


def make_vcf_text(records: Iterable[Tuple[str, int, str, str, str]]) -> str:
    """Return a VCF text containing the standard header and given records.


    records: iterable of (chrom, pos, id, ref, alt)
    """
    lines = [f"{c}\t{p}\t{id_}\t{ref}\t{alt}\t.\tPASS\t." for c, p, id_, ref, alt in records]
    return VCF_HEADER + "\n".join(lines) + "\n"
