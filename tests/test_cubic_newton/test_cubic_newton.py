import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.optimizers.CubicNewton import CubicNewton


@pytest.mark.parametrize(
    "hessian_recalc",
    (None, 1),
)
def test_cubic_newton(hessian_recalc):
    geom = AnaPot.get_geom((-0.25, 1.0, 0))

    opt = CubicNewton(geom, hessian_recalc=hessian_recalc)
    opt.run()

    # calc = geom.calculator
    # calc.plot_opt(opt, show=True)
    assert opt.is_converged
