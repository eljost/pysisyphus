import pytest

from pysisyphus.calculators.Dalton import Dalton
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


@using("dalton")
def test_h2o_opt():
    geom = geom_loader("lib:h2o.xyz", coord_type="redund")
    calc = Dalton(basis="3-21G")
    geom.set_calculator(calc)
    opt = RFOptimizer(geom, thresh="gau_tight")
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == 4
    assert geom.energy == pytest.approx(-75.58595976)
