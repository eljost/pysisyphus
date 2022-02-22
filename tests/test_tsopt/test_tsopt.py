import matplotlib.pyplot as plt
import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.drivers import run_opt
from pysisyphus.helpers import geom_loader
from pysisyphus.tsoptimizers import *
from pysisyphus.testing import using


@pytest.mark.parametrize(
    "opt_cls, ref_cur_cycle",
    [
        pytest.param(TRIM, 9),
        pytest.param(RSIRFOptimizer, 8),
        pytest.param(RSPRFOptimizer, 10),
    ],
)
def test_tshessian_opts(opt_cls, ref_cur_cycle):
    geom = AnaPot.get_geom((-0.6, 2.2, 0.0))

    opt_kwargs = {
        "trust_radius": 0.2,
        "dump": False,
    }
    opt = opt_cls(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == ref_cur_cycle

    # calc = geom.calculator
    # calc.plot_opt(opt, show=True)


@using("pyscf")
@pytest.mark.parametrize(
    "kwargs",
    [
        {"prim_coord": ["BEND", 2, 1, 0]},
        {"hessian_init": "fischer", "rx_coords": [["BEND", 2, 1, 0]]},
    ],
)
def test_tshessianoptimizer_kwargs(kwargs):
    geom = geom_loader("lib:baker_ts/01_hcn.xyz", coord_type="redund")
    geom.set_calculator(PySCF(basis="321g"))

    opt_kwargs = {
        "thresh": "gau_tight",
    }
    opt_kwargs.update(kwargs)
    opt = RSPRFOptimizer(geom, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert geom.energy == pytest.approx(-92.2460426792319)


@using("pyscf")
def test_hessian_ref(this_dir):
    geom = geom_loader("lib:hcn_iso_pm6_near_ts.xyz", coord_type="cart")
    geom.set_calculator(PySCF(basis="321g"))

    hessian_ref = this_dir / "inp_hessian_ref.h5"
    opt = RSPRFOptimizer(geom, hessian_ref=hessian_ref, thresh="gau")
    opt.run()

    assert opt.is_converged
    assert geom.energy == pytest.approx(-92.2460426792319)


def test_iterative_tsopt():
    geom = geom_loader("lib:hcn_iso_pm6_near_ts.xyz", coord_type="cart")
    def calc_getter():
        return PySCF(basis="321g")
    opt_kwargs = {
        "thresh": "gau_loose",
    }
    opt_result = run_opt(geom, calc_getter, "rsprfo", opt_kwargs, iterative=True)
    assert opt_result.opt.is_converged
