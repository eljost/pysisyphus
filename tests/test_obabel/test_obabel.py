import pytest

from pysisyphus.calculators.OBabel import OBabel
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


@using("obabel")
@pytest.mark.parametrize(
    "ff, ref_energy",
    [
        ("mmff94", 0.02585932),
        ("uff", 0.01686154),
        ("gaff", 0.00680588),
        ("ghemical", 0.003055675),
    ],
)
def test_benzene_opt(ff, ref_energy):
    geom = geom_loader("lib:benzene.xyz", coord_type="redund")
    calc = OBabel(ff=ff)
    geom.set_calculator(calc)
    opt = RFOptimizer(geom, thresh="gau_tight")
    opt.run()
    assert geom.energy == pytest.approx(ref_energy)
