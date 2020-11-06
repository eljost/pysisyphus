import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.calculators import DFTBp
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


@using("dftbp")
def test_dftbp_opt():
    geom = geom_loader("lib:h2o.xyz", coord_type="redund")
    calc = DFTBp(parameter="mio-ext")
    geom.set_calculator(calc)

    opt = RFOptimizer(geom, thresh="gau_tight")
    opt.run()

    assert opt.is_converged
    assert geom.energy == pytest.approx(-4.07793793)

