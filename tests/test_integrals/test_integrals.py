import json

import numpy as np
import pytest

from pysisyphus.config import WF_LIB_DIR
from pysisyphus.wavefunction import get_l, Shells
from pysisyphus.wavefunction.cart2sph import cart2sph_coeffs_for


@pytest.mark.parametrize(
    "l_inp, l_ref",
    (
        ("s", 0),
        ("p", 1),
        ("h", 5),
        (17, 17),
        pytest.param(-17, None, marks=pytest.mark.xfail),
    ),
)
def test_get_l(l_inp, l_ref):
    assert get_l(l_inp) == l_ref


@pytest.fixture
def ch4_shells(this_dir):
    fn = WF_LIB_DIR / "turbomole_ch4_def2svp_aomix.in"
    shells = Shells.from_aomix(fn)
    return shells


def test_aomix_shells(ch4_shells):
    assert len(ch4_shells) == 18


@pytest.mark.parametrize("l", range(0, 4))
def test_cart2sph(l):
    Cd = cart2sph_coeffs_for(l, real=True)
    cart_num = (l + 2) * (l + 1) // 2
    sph_num = 2 * l + 1
    assert Cd.shape == (sph_num, cart_num)


@pytest.mark.parametrize(
    "fn",
    (
        "orca_ch4_def2svp.json",
        "orca_ch4_sto3g.json",
        "orca_ch4_tzvpp.json",
        "orca_ch4_qzvpp.json",
    ),
)
def test_orca_spherical_overlaps(fn, this_dir):
    ch4_shells = Shells.from_orca_json(WF_LIB_DIR / fn)
    with open(WF_LIB_DIR / fn) as handle:
        ref_data = json.load(handle)
    S_ref = np.array(ref_data["Molecule"]["S-Matrix"])
    S_sph = ch4_shells.S_sph
    np.testing.assert_allclose(S_sph, S_ref, atol=1e-10)
    mos = ref_data["Molecule"]["MolecularOrbitals"]["MOs"]
    # Also check MO normalization
    C = np.array([mo["MOCoefficients"] for mo in mos]).T
    eye = C.T @ S_sph @ C
    diag = np.diag(eye)
    diag_ref = np.ones_like(diag)
    np.testing.assert_allclose(diag, diag_ref)


def test_cross(this_dir):
    sto3g = Shells.from_orca_json(WF_LIB_DIR / "orca_ch4_sto3g.json")
    tzvpp = Shells.from_orca_json(WF_LIB_DIR / "orca_ch4_tzvpp.json")
    S_sph = sto3g.get_S_sph(tzvpp)
    assert S_sph.shape == (9, 87)


@pytest.mark.skip
def test_orca_dipole_integrals():
    fn = WF_LIB_DIR / "orca_ch4_sto3g.json"
    shells = Shells.from_orca_json(fn)
    C = np.array((0., 0., 0.))
    shells.get_dipole_ints_cart(C)
