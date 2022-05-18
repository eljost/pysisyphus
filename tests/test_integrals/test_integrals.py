import json

import numpy as np
import pytest

from pysisyphus.integrals import get_l, Shells
from pysisyphus.integrals.cart2sph import cart2sph_coeffs_for


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
    fn = this_dir / "methane_def2svp_aomix.in"
    shells = Shells.from_aomix(fn)
    return shells


def test_aomix_shells(ch4_shells):
    assert len(ch4_shells) == 18


@pytest.mark.skip
def test_overlap_matrix(this_dir):
    fn = this_dir / "aomix.in"
    ch4_shells = Shells.from_aomix(fn)
    S_cart = ch4_shells.S_cart
    assert S_cart.shape == (35, 35)


@pytest.mark.parametrize("l", range(0, 4))
def test_cart2sph(l):
    Cd = cart2sph_coeffs_for(l, real=True)
    cart_num = (l + 2) * (l + 1) // 2
    sph_num = 2 * l + 1
    assert Cd.shape == (sph_num, cart_num)


@pytest.mark.parametrize(
    "fn", ("00_ch4.json", "01_ch4_sto3g.json", "02_ch4_tzvpp.json", "04_qzvpp.json")
)
def test_orca_spherical_overlaps(fn, this_dir):
    ch4_shells = Shells.from_orca_json(this_dir / fn)
    with open(this_dir / fn) as handle:
        ref_data = json.load(handle)
    S_ref = np.array(ref_data["Molecule"]["S-Matrix"])
    S_sph = ch4_shells.S_sph
    np.testing.assert_allclose(S_sph, S_ref, atol=1e-10)
