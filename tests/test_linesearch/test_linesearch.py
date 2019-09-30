#!/usr/bin/env python3


import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


def test_linesearch():
    geom = AnaPot.get_geom((0.667, 1.609, 0.))

    rfo_kwargs = {
        # "trust_radius": 0.75,
        "trust_radius": 0.75,
    }
    rfo = RFOptimizer(geom, **rfo_kwargs)
    rfo.run()


    coords = np.array(rfo.coords)
    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    ax.plot(*coords.T[:2], "o-")
    plt.show()


if __name__ == "__main__":
    test_linesearch()
