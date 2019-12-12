#!/usr/bin/env python3

import itertools as it

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


def test_linesearch():
    geom = AnaPot.get_geom((0.667, 1.609, 0.))

    rfo_kwargs = {
        # "trust_radius": 0.75,
        "trust_radius": 0.75,
        "thresh": "gau_tight",
        # "max_cycles": 14,
        # "hessian_recalc": 3,
        # "hessian_init": "calc",
        "line_search": True,
    }
    rfo = RFOptimizer(geom, **rfo_kwargs)
    rfo.run()


    coords = np.array(rfo.coords)
    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    ax.plot(*coords.T[:2], "o-")
    for i, coords in enumerate(coords):
        ax.annotate(f"{i}", coords[:2])
    plt.show()


def run_opt(line_search=False, hessian_recalc=None):
    geom = AnaPot.get_geom((0.667, 1.609, 0.))

    rfo_kwargs = {
        "trust_radius": 0.75,
        "thresh": "gau_tight",
        "hessian_recalc": hessian_recalc,
        "line_search": line_search,
        # "hessian_init": "calc",
    }
    rfo = RFOptimizer(geom, **rfo_kwargs)
    rfo.run()
    conv = rfo.is_converged
    cycs = rfo.cur_cycle
    return (conv, cycs)


def test_line_search():
    ls = (True, False)
    recalc = (None, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    results = list()
    for key in it.product(ls, recalc):
        line_search, hessian_recalc  = key
        conv, cycs = run_opt(line_search, hessian_recalc)
        res = (line_search, hessian_recalc, conv, cycs)
        results.append(res)

    columns = "line_search hessian_recalc converged cycles".split()
    df = pd.DataFrame(results, columns=columns)
    import pdb; pdb.set_trace()

    pass


if __name__ == "__main__":
    test_linesearch()
    # test_line_search()
