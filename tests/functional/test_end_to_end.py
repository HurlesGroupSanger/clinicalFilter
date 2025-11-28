import pytest
import subprocess
import tempfile
import atexit
import os
from tests.helpers.ped_factory import add_absolute_paths


@pytest.mark.functional
def test_pipeline_basic(template_trio_ped, template_ddg2p):

    # Create a temp PED file with absolute path to test VCFs
    tmp_ped = tempfile.NamedTemporaryFile(delete=False, suffix=".ped")
    tmp_ped_path = tmp_ped.name
    tmp_ped.close()
    atexit.register(lambda: os.remove(tmp_ped_path))

    base_path = os.getcwd()
    add_absolute_paths(template_trio_ped, base_path, tmp_ped_path)

    # Run clinical filtering
    cmd = [
        "python",
        "runclinicalfiltering.py",
        "--ped",
        tmp_ped_path,
        "--known-genes",
        template_ddg2p,
        "--outdir",
        f"{base_path}/tests/output",
    ]
    subprocess.check_call(cmd)
