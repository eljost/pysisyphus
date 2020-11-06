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

from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.NCOptimizer import NCOptimizer
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.RSA import RSA


@pytest.mark.parametrize(
    "calc_cls, start, ref_cycle, ref_coords", [
    (AnaPot, (0.667, 1.609, 0.), 12, (1.941, 3.8543, 0.)),
    (AnaPot3, (-0.36, 0.93, 0.), 10, (-1., 0., 0.)),
    (AnaPot4, (-0.50, 3.32, 0.),  12, (-2.2102, 0.3297, 0.)),
    (AnaPotCBM, (-0.32, 0.71, 0.), 10, (-1., 0. ,0.)),
    (CerjanMiller, (-0.46, 1.48, 0.), 10, (0., 0., 0.)),
    (FourWellAnaPot, (1.45, 0.04, 0.), 12, (1.1241, -1.4853, 0.)),
    (LEPSBase, (1.31, 0.82, 0.), 27, (0.74200064, 7.37731942, 0.)),
    (MullerBrownPot, (-0.69, 0.55, 0.), 12, (-0.05, 0.4667, 0.)),
    (Rosenbrock, (-1.00, 1.00, 0.), 40, (1., 1., 0.)),
    ]
)
def test_rfoptimizer(calc_cls, start, ref_cycle, ref_coords):
    geom = calc_cls.get_geom(start)

    print("@Using", calc_cls)

    opt_kwargs = {
        "thresh": "gau_tight",
        "dump": True,
        "overachieve_factor": 2.,
        # "max_micro_cycles": 1,
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
    "opt_cls, opt_kwargs_, ref_cycle", [
    (RFOptimizer, {}, 12),
    (NCOptimizer, {}, 13),
    # LBFGS converged to the saddle point, as the 'hessian' has the
    # wrong eigenvalue structure. Ok, ok we don't have a hessian but
    # you get the idea :)
    pytest.param(LBFGS, {"double_damp": True, "gamma_mult": True, }, 19),
    pytest.param(LBFGS, {"double_damp": True, "gamma_mult": False, }, 19),
    pytest.param(LBFGS, {"double_damp": False}, 19, marks=pytest.mark.xfail),
    (RSA, {}, 17),
    ]
)
def test_optimizers(opt_cls, opt_kwargs_, ref_cycle):
    geom = AnaPot.get_geom((0.667, 1.609, 0.))
    ref_coords = np.array((1.941, 3.8543, 0.))

    opt_kwargs = {
        "thresh": "gau_tight",
        "dump": False,
        "overachieve_factor": 2.,
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
