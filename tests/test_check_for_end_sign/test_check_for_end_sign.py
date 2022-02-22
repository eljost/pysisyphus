from pathlib import Path
import tempfile

import pytest

from pysisyphus.helpers import check_for_end_sign


@pytest.mark.parametrize("check_user", (True, False))
@pytest.mark.parametrize(
    "fn, result", (
    ("stop", True),
    ("converged", True),
    pytest.param("exit", False, marks=pytest.mark.xfail),
    ("nope", False),
))
def test_check_for_end_sign(check_user, fn, result):
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = Path(tmp_dir)
        # touch file
        with open(tmp_dir / fn, "w") as handle:
            pass
        sign_found = check_for_end_sign(check_user=check_user, cwd=tmp_dir)
    assert bool(sign_found) == result
