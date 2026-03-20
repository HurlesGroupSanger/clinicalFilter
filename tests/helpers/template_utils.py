# tests/helpers/template_utils.py
from pathlib import Path
import gzip
from typing import Iterable, Tuple


def copy_template(src: str, dest: Path) -> Path:
    dest.write_bytes(Path(src).read_bytes())
    return dest


def append_vcf_records(dest: Path, records: Iterable[Tuple[str, int, str, str, str]]):
    with dest.open("a", newline="\n") as fh:
        for chrom, pos, id_, ref, alt in records:
            fh.write(f"{chrom}\t{pos}\t{id_}\t{ref}\t{alt}\t.\tPASS\t.\n")


def gzip_file(src: Path, dst: Path) -> Path:
    with src.open("rb") as f_in, gzip.open(dst, "wb") as f_out:
        f_out.writelines(f_in)
    return dst
