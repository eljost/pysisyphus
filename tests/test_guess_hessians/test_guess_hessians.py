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
@pytest.mark.skip_ci
@pytest.mark.parametrize("thresh", ("gau_loose", "gau_tight"))
@pytest.mark.parametrize(
    "hessian_init, ref_cycle_tight, ref_cycle_loose",
    [
        ("calc", 11, 8),
        # Converges to wrong minimum
        # ("unit", 9),
        ("fischer", 16, 13),
        ("lindh", 24, 21),
        ("simple", 22, 4),
        ("swart", 15, 4),
        pytest.param("xtb", 21, 15, marks=[using("pyscf"), using("xtb")]),
        pytest.param("xtb1", 19, 15, marks=[using("pyscf"), using("xtb")]),
    ],
)
def test_guess_hessians(thresh, hessian_init, ref_cycle_tight, ref_cycle_loose):
    """Again, this test seems kind of flakey, esp. the combination of
    XTB and PySCF."""
    geom = geom_loader("lib:h2o2_hf_321g_opt.xyz", coord_type="redund")
    geom.set_calculator(PySCF(basis="def2svp", pal=1))

    print("@\tguess_hessian:", hessian_init)
    opt_kwargs = {
        "hessian_init": hessian_init,
        "thresh": thresh,
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    print("@\tcur_cycle:", opt.cur_cycle)

    assert opt.is_converged
    # ref_cycle = ref_cycle_loose if thresh == "gau_loose" else ref_cycle_tight
    # assert opt.cur_cycle == ref_cycle
    # assert geom.energy == pytest.approx(-150.65298169)  # tight
    assert geom.energy == pytest.approx(-150.652, abs=1e-3)  # loose


@using("pyscf")
@pytest.mark.parametrize(
    "hessian_init",
    (
        "calc",
        "unit",
        "fischer",
        "lindh",
        "simple",
        "swart",
        pytest.param("xtb", marks=[using("pyscf"), using("xtb")]),
        pytest.param("xtb1", marks=[using("pyscf"), using("xtb")]),
    ),
)
def test_gen_guess_hessians(hessian_init):
    geom = geom_loader("lib:h2o2_hf_321g_opt.xyz", coord_type="redund")
    geom.set_calculator(PySCF(basis="def2svp", pal=1))
    opt = RFOptimizer(geom, hessian_init=hessian_init)
    opt.prepare_opt()
    assert opt.H is not None


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
