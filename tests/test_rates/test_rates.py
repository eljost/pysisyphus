import pytest

from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.drivers import eyring_rate


def test_eyring():
    barrier = 104.5 / AU2KJPERMOL  # about 25 kcal mol⁻¹
    T = 298.15
    rate = eyring_rate(barrier, temperature=T)
    assert rate == pytest.approx(3.059526e-06 )
