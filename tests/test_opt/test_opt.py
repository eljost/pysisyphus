import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.AnaPot3 import AnaPot3
from pysisyphus.calculators.AnaPot4 import AnaPot4
from pysisyphus.calculators.AnaPotCBM import AnaPotCBM
from pysisyphus.calculators.CerjanMiller import CerjanMiller
from pysisyphus.calculators.FourWellAnaPot import FourWellAnaPot
from pysisyphus.calculators.LEPSBase import LEPSBase
from pysisyphus.calculators.MullerBrownSympyPot import MullerBrownPot
from pysisyphus.calculators.Rosenbrock import Rosenbrock

from pysisyphus.optimizers.BFGS import BFGS
from pysisyphus.optimizers.closures import lbfgs_closure, modified_broyden_closure
from pysisyphus.optimizers.ConjugateGradient import ConjugateGradient
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.NCOptimizer import NCOptimizer
from pysisyphus.optimizers.restrict_step import scale_by_max_step
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.RSA import RSA
from pysisyphus.optimizers.StabilizedQNMethod import StabilizedQNMethod


@pytest.mark.parametrize(
    "calc_cls, start, ref_cycle, ref_coords",
    [
        (AnaPot, (0.667, 1.609, 0.0), 18, (1.941, 3.8543, 0.0)),
        (AnaPot3, (-0.36, 0.93, 0.0), 10, (-1.0, 0.0, 0.0)),
        (AnaPot4, (-0.50, 3.32, 0.0), 11, (-2.2102, 0.3297, 0.0)),
        (AnaPotCBM, (-0.32, 0.71, 0.0), 10, (-1.0, 0.0, 0.0)),
        (CerjanMiller, (-0.46, 1.48, 0.0), 10, (0.0, 0.0, 0.0)),
        (FourWellAnaPot, (1.45, 0.04, 0.0), 10, (1.1241, -1.4853, 0.0)),
        (LEPSBase, (1.31, 0.82, 0.0), 27, (0.74200064, 7.17588688, 0.0)),
        (MullerBrownPot, (-0.69, 0.55, 0.0), 12, (-0.05, 0.4667, 0.0)),
        (Rosenbrock, (-1.00, 1.00, 0.0), 43, (1.0, 1.0, 0.0)),
    ],
)
def test_rfoptimizer(calc_cls, start, ref_cycle, ref_coords):
    geom = calc_cls.get_geom(start)

    print("@Using", calc_cls)

    opt_kwargs = {
        "thresh": "gau_tight",
        "dump": True,
        "overachieve_factor": 2.0,
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    # geom.calculator.plot_opt(opt, show=True)
    # import matplotlib.pyplot as plt
    # calc = geom.calculator
    # calc.plot()
    # coords = np.array(opt.coords)
    # ax = calc.ax
    # ax.plot(*coords.T[:2], "ro-")
    # plt.show()

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle

    ref_coords = np.array(ref_coords)
    diff = ref_coords - geom.coords
    diff_norm = np.linalg.norm(diff)
    print(f"@\tnorm(diff)={diff_norm:.8f}")
    assert diff_norm < 6e-5

    print("@\tFinal coords", geom.coords)


@pytest.mark.parametrize(
    "opt_cls, opt_kwargs_, ref_cycle",
    [
        (RFOptimizer, {}, 18),
        (RFOptimizer, {"adapt_step_func": True}, 18),
        (NCOptimizer, {}, 13),
        # LBFGS converged to the saddle point, as the 'hessian' has the
        # wrong eigenvalue structure. Ok, ok we don't have a hessian but
        # you get the idea :)
        pytest.param(
            LBFGS,
            {
                "double_damp": True,
                "gamma_mult": True,
            },
            19,
        ),
        pytest.param(
            LBFGS,
            {
                "double_damp": True,
                "gamma_mult": False,
            },
            19,
        ),
        pytest.param(LBFGS, {"double_damp": False}, 19, marks=pytest.mark.xfail),
        pytest.param(BFGS, {}, 10, marks=pytest.mark.xfail),
        (RSA, {}, 17),
        pytest.param(ConjugateGradient, {"formula": "PR"}, 170),
        # Something goes terribly wrong here ...
        pytest.param(StabilizedQNMethod, {"bio": False}, 10, marks=pytest.mark.xfail),
    ],
)
def test_optimizers(opt_cls, opt_kwargs_, ref_cycle):
    geom = AnaPot.get_geom((0.667, 1.609, 0.0))
    ref_coords = np.array((1.941, 3.8543, 0.0))

    opt_kwargs = {
        "thresh": "gau_tight",
        "dump": False,
        "overachieve_factor": 2.0,
        # Add one, so max_cycles is always bigger than ref_cycle
        "max_cycles": max(ref_cycle + 1, 50),
    }
    opt_kwargs.update(opt_kwargs_)
    opt = opt_cls(geom, **opt_kwargs)
    opt.run()

    # import matplotlib.pyplot as plt
    # calc = geom.calculator
    # calc.plot()
    # coords = np.array(opt.coords)
    # ax = calc.ax
    # ax.plot(*coords.T[:2], "ro-")
    # plt.show()

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle

    diff = ref_coords - geom.coords
    diff_norm = np.linalg.norm(diff)
    print(f"@\tnorm(diff)={diff_norm:.8f}")
    assert diff_norm < 6e-5

    # print("@\tFinal coords", geom.coords)


def test_thresh_never():
    geom = AnaPot.get_geom((0.667, 1.609, 0.0))

    opt_kwargs = {
        "thresh": "never",
        "assert_min_step": False,
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    # Restrict max_cycles; will crash otherwise, as the BFGS upadte becomes
    # faulty for vanishing step and gradient differences.
    opt.max_cycles = 20
    opt.run()

    assert geom.energy == pytest.approx(0.98555442)
    norm = np.linalg.norm(geom.forces)
    assert norm == pytest.approx(0.0)


@pytest.mark.parametrize(
    "closure_func",
    (
        lbfgs_closure,
        # pytest.param(modified_broyden_closure, marks=pytest.mark.skip),
        modified_broyden_closure,
    ),
)
def test_opt_closure(closure_func):
    ref_coords = np.array((-1.05280566, 1.02776684, 0.0))
    geom = AnaPot.get_geom((-0.8, 1.4, 0.0))

    def forces_getter(coords):
        return geom.get_energy_and_cart_forces_at(coords)["forces"]

    def restrict_step(coords, step):
        return scale_by_max_step(step, 0.1)

    step_func = closure_func(forces_getter, restrict_step=restrict_step)
    for _ in range(50):
        step, forces = step_func(geom.coords)
        force_norm = np.linalg.norm(forces)
        print(_, force_norm)
        if force_norm <= 1e-6:
            print("Converged")
            break
        geom.coords = geom.coords + step
    else:
        raise AssertionError
    assert np.linalg.norm(geom.coords - ref_coords) == pytest.approx(0.0, abs=1e-4)
