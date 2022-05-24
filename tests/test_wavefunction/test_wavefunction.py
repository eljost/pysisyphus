import numpy as np
import pytest

from pysisyphus.config import WF_LIB_DIR
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction.pop_analysis import (
    mulliken_charges_from_wf,
    iao_charges_from_wf,
)


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
    charges = mulliken_charges_from_wf(wf)
    np.testing.assert_allclose(
        charges,
        (-0.27748112, 0.06937027, 0.06937027, 0.06937027, 0.06937027),
        atol=2e-5,
    )


@pytest.mark.parametrize(
    "fn",
    (
        "iao_ref.json",
        "iao_ref_uhf.json",
    ),
)
def test_orca_iao(fn):
    wf = Wavefunction.from_orca_json(WF_LIB_DIR / fn)
    charges = iao_charges_from_wf(wf)
    np.testing.assert_allclose(
        charges,
        (-0.523749, 0.130937, 0.130937, 0.130937, 0.130937),
        atol=2e-5,
    )


@pytest.mark.parametrize(
    "fn, ref_dip_mom",
    (
        ("orca_dipmom.json", (0.65100, -0.15774, -0.67608)),
        ("orca_dipmom_uhf.json", (0.00000, 0.00000, -0.34612)),
    ),
)
def test_orca_dipole_moments(fn, ref_dip_mom):
    wf = Wavefunction.from_orca_json(WF_LIB_DIR / fn)
    dip_mom = wf.dipole_moment()
    np.testing.assert_allclose(dip_mom, ref_dip_mom, atol=1e-5)
