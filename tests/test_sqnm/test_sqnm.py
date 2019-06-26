#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_library
from pysisyphus.init_logging import init_logging
from pysisyphus.optimizers.StabilizedQNMethod import StabilizedQNMethod
from pysisyphus.optimizers.StabilizedQNMethod2 import StabilizedQNMethod as SQNM2
from pysisyphus.optimizers.SQNM_ref import StabilizedQNMethod as SQNM3
from pysisyphus.calculators.XTB import XTB


def test_sqnm():
    # geom = AnaPot.get_geom((0, 4, 0))
    geom = AnaPot.get_geom((-0.8, 1.73, 0))

    opt_kwargs = {
        "max_cycles": 15,
        # "max_cycles": 5,
        "eps": 1e-4,
        "E_thresh": 1e-4,
        "alpha": 0.5,
        "hist_max": 10,
        "thresh": "gau_tight",
        "trust_radius": 0.1,
    }
    opt = StabilizedQNMethod(geom, **opt_kwargs)
    # opt = SQNM3(geom, **opt_kwargs)
    opt.run()
    c = np.array(opt.coords)
    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    ax.plot(c[:,0], c[:,1], "o-")

    plt.show()


def test_sqnm_bio_mode():
    atoms = "h h".split()
    coords = ((0, 0, 0), (0, 0, 1))
    geom = Geometry(atoms, coords)
    xtb = XTB()
    geom.set_calculator(xtb)

    opt = StabilizedQNMethod(geom)
    cur_grad = geom.gradient
    stretch_grad, rem_grad = opt.bio_mode(cur_grad)
    # In H2 there is only one bond, so stretch_gradient == cur_grad and the
    # remainder should be the zero vector.
    np.testing.assert_allclose(stretch_grad, cur_grad)
    np.testing.assert_allclose(np.linalg.norm(rem_grad), 0.)
    return


def test_sqnm_xtb():
    geom = geom_from_library("split.image_021.xyz")
    # geom = geom_from_library("split.image_021.xyz", coord_type="redund")
    xtb = XTB()
    geom.set_calculator(xtb)

    opt_kwargs = {
        "max_cycles": 93,
        "eps": 1e-4,
        "E_thresh": 1e-6,
        "alpha": 0.5,
        "alpha_stretch": 0.5,
        "hist_max": 10,
        "dump": True,
        "bio": False,
        "trust_radius": 0.1,
    }
    opt = StabilizedQNMethod(geom, **opt_kwargs)
    # # opt = SQNM2(geom, **opt_kwargs)
    # opt = SQNM3(geom, **opt_kwargs)

    # from pysisyphus.optimizers.RFOptimizer import RFOptimizer
    # opt = RFOptimizer(geom)

    # from pysisyphus.optimizers.SteepestDescent import SteepestDescent
    # opt = SteepestDescent(geom, max_cycles=150)

    # from pysisyphus.optimizers.BFGS import BFGS
    # opt = BFGS(geom, max_cycles=150)

    opt.run()


if __name__ == "__main__":
    # test_sqnm()
    # test_sqnm_bio_mode()
    test_sqnm_xtb()
