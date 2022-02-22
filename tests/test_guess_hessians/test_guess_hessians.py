import pytest

from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_loader, do_final_hessian
from pysisyphus.optimizers.guess_hessians import ts_hessian
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer
from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer
from pysisyphus.testing import using

import numpy as np


@using("pyscf")
@pytest.mark.parametrize(
    "hessian_init, ref_cycle",
    [
        ("calc", 11),
        # Converges to wrong minimum
        # ("unit", 9),
        ("fischer", 16),
        ("lindh", 24),
        ("simple", 22),
        ("swart", 15),
        pytest.param("xtb", 21, marks=[using("pyscf"), using("xtb")]),
        pytest.param("xtb1", 19, marks=[using("pyscf"), using("xtb")]),
    ],
)
def test_guess_hessians(hessian_init, ref_cycle):
    geom = geom_loader("lib:h2o2_hf_321g_opt.xyz", coord_type="redund")
    geom.set_calculator(PySCF(basis="def2svp", pal=2))

    print("@\tguess_hessian:", hessian_init)
    opt_kwargs = {
        "hessian_init": hessian_init,
        "thresh": "gau_tight",
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    print("@\tcur_cycle:", opt.cur_cycle)

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle
    assert geom.energy == pytest.approx(-150.65298169)


@using("pyscf")
def test_ts_hessian():
    H = np.diag((1, 0.5, 0.25))
    tsh = ts_hessian(H, (0, 1))
    w, v = np.linalg.eigh(tsh)
    np.testing.assert_allclose(w, (-0.445194, 0.070194, 0.25), atol=1e-7)


@using("pyscf")
@pytest.mark.parametrize(
    "tsopt_cls, ref_cycle",
    [
        (RSPRFOptimizer, 10),
        # (RSIRFOptimizer, 12),
    ],
)
def test_ts_hessian_opt(tsopt_cls, ref_cycle):
    geom = geom_loader("lib:baker_ts/01_hcn.xyz", coord_type="redund")
    geom.set_calculator(PySCF(basis="321g"))

    opt_kwargs = {
        "hessian_init": "fischer",
        "dump": True,
        "rx_coords": (("BEND", 2, 1, 0),),
        "thresh": "gau_tight",
        "trust_max": 0.3,
    }
    opt = tsopt_cls(geom, **opt_kwargs)
    opt.run()
    do_final_hessian(geom, write_imag_modes=True)

    assert opt.is_converged
    # 11 without linesearch
    assert opt.cur_cycle == ref_cycle
    assert geom.energy == pytest.approx(-92.2460426792319)
