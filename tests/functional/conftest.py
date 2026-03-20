# tests/functional/conftest.py
import pytest
from pathlib import Path

# from tests.helpers.vcf_factory import make_vcf_text
# from tests.helpers.ped_factory import make_ped_text
# from tests.helpers.template_utils import copy_template, append_vcf_records, gzip_file


TEMPLATES = Path(__file__).parent.parent / "data" / "templates"


@pytest.fixture
def template_trio_ped():
    return str(TEMPLATES / "trio.ped")


@pytest.fixture
def template_ddg2p():
    return str(TEMPLATES / "ddg2p.tsv")


# @pytest.fixture
# def make_case_from_template(tmp_path):
#     def _make(template_vcf, template_ped, *, records_to_add=None, gzip_output=False, ped_members_override=None):
#         vcf_out = tmp_path / "case.vcf"
#         ped_out = tmp_path / "case.ped"
#         copy_template(template_vcf, vcf_out)
#         copy_template(template_ped, ped_out)

#         if records_to_add:
#             append_vcf_records(vcf_out, records_to_add)

#         if ped_members_override:
#             ped_out.write_text(make_ped_text(ped_members_override["family_id"], ped_members_override["members"]))

#         if gzip_output:
#             gz = tmp_path / "case.vcf.gz"
#             gzip_file(vcf_out, gz)
#             return {"vcf": str(gz), "ped": str(ped_out)}

#         return {"vcf": str(vcf_out), "ped": str(ped_out)}

#     return _make
