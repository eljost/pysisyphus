#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.helpers import geom_from_library
from pysisyphus.init_logging import init_logging
from pysisyphus.optimizers.StabilizedQNMethod import StabilizedQNMethod
from pysisyphus.calculators.XTB import XTB


def test_sqnm():
    # geom = AnaPot.get_geom((0, 4, 0))
    geom = AnaPot.get_geom((-0.8, 1.73, 0))

    opt_kwargs = {
        "max_cycles": 15,
        # "max_cycles": 5,
        "eps": 2e-4,
        "E_thresh": 1e-4,
        "alpha": 0.01,
        "hist_max": 5,
    }
    opt = StabilizedQNMethod(geom, **opt_kwargs)
    opt.run()
    c = np.array(opt.coords)
    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    ax.plot(c[:,0], c[:,1], "o-")

    plt.show()


def test_sqnm_xtb():
    geom = geom_from_library("split.image_021.xyz", coord_type="redund")
    xtb = XTB()
    geom.set_calculator(xtb)

    opt_kwargs = {
        "max_cycles": 150,
        "eps": 1e-4,
        "E_thresh": 1e-6,
        "alpha": 0.5,
        "hist_max": 10,
        "dump": True,
    }
    opt = StabilizedQNMethod(geom, **opt_kwargs)
    opt.bio_mode()
    return

    # from pysisyphus.optimizers.RFOptimizer import RFOptimizer
    # opt = RFOptimizer(geom)

    # from pysisyphus.optimizers.SteepestDescent import SteepestDescent
    # opt = SteepestDescent(geom, max_cycles=150)

    # from pysisyphus.optimizers.BFGS import BFGS
    # opt = BFGS(geom, max_cycles=150)

    opt.run()


if __name__ == "__main__":
    # test_sqnm()
    test_sqnm_xtb()
