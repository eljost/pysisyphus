import pytest

from pysisyphus.calculators.Rastrigin import Rastrigin
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


def test_rastrigin_minimum():
    geom = Rastrigin().get_minima()[0]

    # calc = geom.calculator
    # calc.plot(show=True)

    opt = RFOptimizer(geom, thresh="gau_tight")
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == 0
    assert geom.energy == pytest.approx(0.)
