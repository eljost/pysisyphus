import pytest

from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.drivers import eyring_rate, harmonic_tst_rate


T = 298.15
BARRIER = 104.5 / AU2KJPERMOL # about 25 kcal mol⁻¹


def test_eyring():
    rate = eyring_rate(BARRIER, temperature=T)
    assert rate == pytest.approx(3.059526e-06 )


@pytest.mark.parametrize(
    "rs_part_func, ts_part_func, ref_rate", (
        (1., 1., 3.059526e-6),
    )
)
def test_harmonic_tst(rs_part_func, ts_part_func, ref_rate):
    rate = harmonic_tst_rate(BARRIER, T, rs_part_func, ts_part_func)
    assert rate == pytest.approx(ref_rate)
