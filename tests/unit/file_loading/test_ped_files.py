import pytest
from file_loading.ped_files import openped


TEMPLATES_DIR = "tests/data/templates"


@pytest.mark.parametrize(
    "ped_file,proband_list",
    [(f"{TEMPLATES_DIR}/trio.ped", None), (f"{TEMPLATES_DIR}/singleton.ped", None)],
)
def test_openped(ped_file, proband_list):

    families = openped(ped_file, proband_list)
    assert len(families) == 1
