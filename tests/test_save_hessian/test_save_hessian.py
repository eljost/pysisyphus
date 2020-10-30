import pytest

from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.irc import EulerPC
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using


@using("pyscf")
def test_opt_load_save_hessian(this_dir):
    geom = geom_loader("lib:acetaldehyd.xyz", coord_type="redund")
    calc = PySCF(basis="321g", pal=2)
    geom.set_calculator(calc)

    opt_kwargs = {
        "hessian_recalc": 1,
        "hessian_init": this_dir / "inp_hess_calc_cyc_0.h5",
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == 2
    assert geom.energy == pytest.approx(-152.05524620313963)


@using("pyscf")
def test_irc_load_save_hessian(this_dir):
    geom = geom_loader("lib:hcn_iso_hf_sto3g_ts_opt.xyz")
    calc = PySCF(basis="sto3g", pal=2)
    geom.set_calculator(calc)

    irc_kwargs = {
        "rms_grad_thresh": 1e-3,
        "hessian_init": this_dir / "inp_hess_init_irc.h5",
        "hessian_recalc": 2,
        "max_cycles": 5,
    }
    irc = EulerPC(geom, **irc_kwargs)
    irc.run()

    assert irc.forward_cycle == 4
    assert irc.backward_cycle == 4
