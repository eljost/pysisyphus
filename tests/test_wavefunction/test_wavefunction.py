import numpy as np
import pytest

from pysisyphus.config import WF_LIB_DIR
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction.pop_analysis import mulliken_from_wf


@pytest.mark.parametrize(
    "fn, unrestricted",
    (
        ("orca_ch4_sto3g.json", False),
        ("orca_ch4_sto3g_uhf.json", True),
    ),
)
def test_orca_json(fn, unrestricted):
    wf = Wavefunction.from_orca_json(WF_LIB_DIR / fn)
    assert wf.unrestricted == unrestricted
    assert wf.occ == (5, 5)


@pytest.mark.parametrize(
    "fn",
    ("orca_ch4_sto3g.json", "orca_ch4_sto3g_uhf.json"),
)
def test_orca_mulliken(fn):
    wf = Wavefunction.from_orca_json(WF_LIB_DIR / fn)
    charges = mulliken_from_wf(wf)
    np.testing.assert_allclose(
        charges,
        (-0.27748112, 0.06937027, 0.06937027, 0.06937027, 0.06937027),
        atol=2e-5,
    )
