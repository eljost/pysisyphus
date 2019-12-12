#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.helpers import geom_from_library
from pysisyphus.calculators.CerjanMiller import CerjanMiller
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer
from pysisyphus.calculators.XTB import XTB
from pysisyphus.init_logging import init_logging


def test_rsprfo_hcn_ts_xtb():
    geom = geom_from_library("hcn_iso_ts.xyz", coord_type="redund")
    xtb = XTB()
    geom.set_calculator(xtb)

    opt_kwargs = {
        "thresh": "gau_tight",
        "max_micro_cycles": 1,
    }
    opt = RSPRFOptimizer(geom, **opt_kwargs)
    opt.run()
    assert opt.is_converged
    assert opt.cur_cycle == 7


def test_prfo_analytical():
    geom = CerjanMiller.get_geom((0.559714, -0.4885, 0))
    opt_kwargs = {
        "trust_max": 0.1,
    }
    # Without RS
    opt = RSPRFOptimizer(geom, max_micro_cycles=1, **opt_kwargs)
    opt.run()
    assert opt.is_converged
    assert opt.cur_cycle == 9

    rs_geom = CerjanMiller.get_geom((0.559714, -0.4885, 0))
    rs_opt = RSPRFOptimizer(rs_geom, **opt_kwargs)
    rs_opt.run()
    assert rs_opt.is_converged
    assert rs_opt.cur_cycle == 8

    # cs = np.array(opt.coords)
    # rs_cs = np.array(rs_opt.coords)
    # calc = geom.calculator
    # calc.plot()
    # ax = calc.ax
    # ax.plot(*cs.T[:2], "ro-", label="no RS")
    # ax.plot(*rs_cs.T[:2], "go-", label="RS")
    # ax.legend()
    # plt.show()


if __name__ == "__main__":
    test_rsprfo_hcn_ts_xtb()
    init_logging()
    test_prfo_analytical()
