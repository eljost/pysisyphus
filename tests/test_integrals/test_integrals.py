import pytest

from pysisyphus.integrals import get_l, Shells


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
    fn = this_dir  / "methane_def2svp_aomix.in"
    shells = Shells.from_aomix(fn)
    return shells


def test_aomix_shells(ch4_shells):
    assert len(ch4_shells) == 18


@pytest.mark.skip
def test_overlap_matrix(this_dir):
    fn = this_dir / "aomix.in"
    ch4_shells = Shells.from_aomix(fn)
    S_cart = ch4_shells.S_cart
    import numpy as np
    S_cart_ref = np.load("/home/johannes/Code/overlaps/orbkit/S_sto3g.npy")
    import pdb; pdb.set_trace()  # fmt: skip
    assert S_cart.shape == (35, 35)
