import pytest

from pysisyphus.calculators import Gaussian16
from pysisyphus.testing import using


@pytest.mark.parametrize("route", ("symmetry", "nosymm", "force", "opt", "freq", "irc"))
def test_invalid_keywords(route):
    with pytest.raises(AssertionError):
        Gaussian16(route)


@using("gaussian16")
def test_all_energies():
    geom = Gaussian16.geom_from_fn("lib:h2o.xyz", route="hf def2svp td=(nstates=2)")
    all_energies = geom.all_energies
    len(all_energies) == 3
    assert all_energies[2] == pytest.approx(-75.5569539)
