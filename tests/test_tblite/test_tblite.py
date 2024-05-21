import pytest

from pysisyphus.calculators.TBLite import TBLite
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using


@using("tblite")
@pytest.mark.parametrize(
    "gfn, ref_energy",
    [
        # GFN0 and GFNFF don't seem to be supported by tblite
        # (0, -4.36679493),
        (1, -5.76865989),
        (2, -5.07043090),
        # ("ff", -0.327181905),
    ],
)
def test_tblite(gfn, ref_energy):
    geom = geom_loader("lib:h2o.xyz")
    calc = TBLite(gfn=gfn)
    geom.set_calculator(calc)

    energy = geom.energy
    forces = geom.forces
    assert energy == pytest.approx(ref_energy)
