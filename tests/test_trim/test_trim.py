#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.tsoptimizers.TRIM import TRIM
from pysisyphus.tsoptimizers.RSIRFO import RSIRFO


def test_trim():
    geom = AnaPot.get_geom((-0.6, 2.2, 0.))

    trim = TRIM(geom, trust_radius=0.2)
    trim.run()

    # assert trim.cur_cycle == 17
    # assert trim.is_converged

    # cs = np.array(trim.coords)
    # calc = geom.calculator
    # calc.plot()
    # ax = calc.ax
    # ax.plot(*cs.T[:2])
    # plt.show()


def test_rsirfo():
    geom = AnaPot.get_geom((-0.6, 2.2, 0.))

    trim = RSIRFO(geom, trust_radius=0.2)
    trim.run()

    # assert trim.cur_cycle == 17
    # assert trim.is_converged

    cs = np.array(trim.coords)
    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    ax.plot(*cs.T[:2], "o-")
    plt.show()


if __name__ == "__main__":
    test_trim()
    test_rsirfo()
