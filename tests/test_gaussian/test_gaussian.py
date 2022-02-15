import pytest

from pysisyphus.calculators import Gaussian16


@pytest.mark.parametrize("route", ("symmetry", "nosymm", "force", "opt", "freq", "irc"))
def test_invalid_keywords(route):
    with pytest.raises(AssertionError):
        Gaussian16(route)
