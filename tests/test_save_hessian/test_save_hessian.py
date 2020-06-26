import pytest

from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


@using("pyscf")
def test_save_hessian(this_dir):
    geom = geom_loader("lib:acetaldehyd.xyz", coord_type="redund")
    calc = PySCF(basis="321g", pal=2)
    geom.set_calculator(calc)

    opt_kwargs = {
        "hessian_recalc": 1,
        "hessian_init": this_dir / "hess_calc_cyc_0.h5",
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == 1
    assert geom.energy == pytest.approx(-152.05524620313963)
