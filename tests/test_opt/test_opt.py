#!/usr/bin/env python3

import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.AnaPot3 import AnaPot3
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


@pytest.mark.parametrize(
    "calc, start, ref_cycle", [
    (AnaPot, (0.667, 1.609, 0.), 11),
    (AnaPot3, (0.1, 0.6, 0.), 11),
    ]
)
def test_rfoptimizer(calc, start, ref_cycle):
    geom = calc.get_geom(start)

    opt_kwargs = {
        "thresh": "gau_tight",
        "dump": False,
    }
    opt = RFOptimizer(geom)
    opt.run()
    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle

    # import matplotlib.pyplot as plt
    # import numpy as np
    # calc = geom.calculator
    # calc.plot()
    # coords = np.array(opt.coords)
    # ax = calc.ax
    # ax.plot(*coords.T[:2], "ro-")
    # plt.show()
