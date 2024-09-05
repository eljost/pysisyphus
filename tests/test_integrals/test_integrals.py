import json

import numpy as np
import pytest

from pysisyphus.config import WF_LIB_DIR
from pysisyphus.io.fchk import parse_fchk
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


@pytest.mark.parametrize(
    "fn",
    (
        "orca_h2o_sto3g.json",
        "orca_h2o_def2svp.json",
    ),
)
def test_orca_one_el_integrals(fn):
    fn = WF_LIB_DIR / fn
    shells = Shells.from_orca_json(fn, backend="python")
    V = shells.V_sph
    T = shells.T_sph
    H = T + V

    """
    with np.printoptions(suppress=True, precision=5, linewidth=140):
        print()
        print("H=T+V")
        print(H)
        print()
        print("T")
        print(T)
        print("V")
        print(V)
    """

    with open(fn) as handle:
        ref_data = json.load(handle)["Molecule"]

    T_ref = ref_data["T-Matrix"]
    np.testing.assert_allclose(T, T_ref, atol=1e-8)

    H_ref = ref_data["H-Matrix"]
    np.testing.assert_allclose(H, H_ref, atol=1e-8)


@pytest.mark.parametrize(
    "fchk",
    ("g16_ch4_qzvpp.fchk",),
)
def test_fchk_overlaps(fchk):
    shells = Shells.from_fchk(WF_LIB_DIR / fchk)
    # Recover reference overlap matrix from MO coefficients
    data = parse_fchk(WF_LIB_DIR / fchk)
    C = np.array(data["Alpha MO coefficients"], dtype=float)
    bf_num = int(np.sqrt(C.size))
    # Transpose, so MOs are in columns
    C = C.reshape(-1, bf_num).T
    C_inv = np.linalg.pinv(C, rcond=1e-8)
    S_ref = C_inv.T @ C_inv
    S = shells.S_sph
    if S.shape != S_ref.shape:
        S = shells.S_cart
    np.testing.assert_allclose(S, S_ref, atol=5e-8)
