#!/usr/bin/env python3


import matplotlib.pyplot as plt
import numpy as np


from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer


def test_pred():
    geom = AnaPot.get_geom((-0.174, 1.27, 0.))
    calc = geom.calculator
    trst = 0.05
    opt_kwargs = {
        "trust_min": trst,
        "trust_radius": trst,
        "trust_max": 0.3,
        "hessian_recalc": 1,
        "max_micro_cycles": 1,
    }
    opt = RSPRFOptimizer(geom, **opt_kwargs)
    opt.run()
    coords = np.array(opt.coords)

    calc.plot()
    ax = calc.ax
    ax.plot(*coords.T[:2], "ro-")
    plt.show()


if __name__ == "__main__":
    test_pred()
