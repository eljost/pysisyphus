#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.tsoptimizers.dimerv2 import dimer_method


def run_dimer():
    x0 = (-0.5767, 1.6810, 0)
    geom0 = AnaPot.get_geom(x0)

    N = np.array((-0.9, 0.43, 0))
    N /= np.linalg.norm(N)
    R = 0.125
    coords = dimer_method(geom0, N, R, AnaPot, max_step=0.6, max_cycles=5)

    coords = np.array(coords)
    pot = AnaPot()
    pot.plot()
    ax = pot.ax
    for i, rot_cycle in enumerate(coords):
        ax.plot(*rot_cycle.T[:2], "o-", label=f"Cycle {i:02d}")
    ax.legend()
    plt.show()

    # Reference
    # results_-7437746119687769865.pickle


if __name__ == "__main__":
    run_dimer()
