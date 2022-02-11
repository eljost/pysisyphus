import numpy as np
import pytest

from pysisyphus.calculators import XTB
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader
from pysisyphus.optimizers.hessian_updates import (
    bfgs_update,
    damped_bfgs_update,
    double_damp,
    sr1_update,
    psb_update,
    flowchart_update,
    mod_flowchart_update,
    bofill_update,
)
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.testing import using
from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer


@pytest.mark.parametrize(
    "update_func",
    [
        bfgs_update,
        damped_bfgs_update,
        flowchart_update,
        mod_flowchart_update,
        bofill_update,
    ],
)
def test_hessian_updates(update_func):
    N = 3
    dx = np.ones(N)
    dg = np.ones(N)
    H = np.arange(N * N).reshape(-1, N)
    dH = update_func(H, dx, dg)


@using("pyscf")
@pytest.mark.parametrize(
    "hessian_update",
    [
        "bfgs",
        "none",
    ],
)
def test_no_hessian_update(hessian_update):
    geom = geom_loader("lib:h2o.xyz")
    calc = PySCF(basis="sto3g", pal=2)
    geom.set_calculator(calc)
    opt = RFOptimizer(geom, thresh="gau", hessian_update=hessian_update)
    opt.run()

    assert geom.energy == pytest.approx(-74.96590119)


@using("xtb")
@pytest.mark.parametrize(
    "hessian_update",
    (
        "bofill",
        "ts_bfgs",
        "ts_bfgs_org",
        "ts_bfgs_rev",
    ),
)
def test_ts_hessian_update(this_dir, hessian_update):
    geom = geom_loader("lib:tsbfgs_init.xyz", coord_type="redund")
    calc = XTB(pal=6)

    geom.set_calculator(calc)

    opt = RSIRFOptimizer(
        geom,
        hessian_init=this_dir / "tsbfgs_init_hess.h5",
        hessian_update=hessian_update,
    )
    opt.run()
    assert opt.is_converged
    assert geom.energy == pytest.approx(-17.81225910)


@using("pyscf")
@pytest.mark.parametrize(
    "hessian_update",
    (
        "ts_bfgs",
        "ts_bfgs_org",
        "ts_bfgs_rev",
    ),
)
def test_ts_hessian_update_hcn(hessian_update):
    geom = geom_loader("lib:hcn_iso_pm6_near_ts.xyz", coord_type="redund")
    calc = PySCF(basis="321g", pal=2)

    geom.set_calculator(calc)

    opt = RSIRFOptimizer(
        geom,
        hessian_update=hessian_update,
    )
    opt.run()
    assert opt.is_converged
    assert geom.energy == pytest.approx(-92.24603904)
