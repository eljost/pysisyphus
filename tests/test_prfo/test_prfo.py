#!/usr/bin/env python3


import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.CerjanMiller import CerjanMiller
from pysisyphus.helpers import geom_from_library
from pysisyphus.tsoptimizers.PRFOptimizer import PRFOptimizer
from pysisyphus.optimizers.RSRFOptimizer import RSRFOptimizer
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer
from pysisyphus.tsoptimizers.PRFOptimizer import PRFOptimizer
from pysisyphus.calculators.XTB import XTB


def check_eigvals(H):
    w, v = np.linalg.eigh(H)
    neg_inds = w < -1e-8
    neg_num = neg_inds.sum()
    eigval_str = np.array2string(w[neg_inds], precision=6)
    print(f"Found {neg_num} negative eigenvalue(s): {eigval_str}")


def test_rsprfo_hcn_ts_xtb():
    geom = geom_from_library("hcn_iso_ts.xyz", coord_type="redund")
    xtb = XTB()
    geom.set_calculator(xtb)

    print("Start")
    check_eigvals(geom.hessian)

    opt_kwargs = {
        "thresh": "gau_tight",
    }
    opt = PRFOptimizer(geom, **opt_kwargs)
    opt.run()
    assert opt.is_converged
    assert opt.cur_cycle == 7

    print()
    print("End")
    check_eigvals(geom.hessian)


def test_prfo_analytical():
    geom = CerjanMiller.get_geom((0.559714, -0.4885, 0))
    opt_kwargs = {
        "trust_max": 0.1,
        # "hessian_recalc": 1,
    }
    opt = PRFOptimizer(geom, **opt_kwargs)
    opt.run()
    assert opt.is_converged
    assert opt.cur_cycle == 9
    # cs = np.array(opt.coords)
    # calc = geom.calculator
    # calc.plot()
    # ax = calc.ax
    # ax.plot(*cs.T[:2], "ro-")
    # plt.show()


if __name__ == "__main__":
    test_rsprfo_hcn_ts_xtb()
    test_prfo_analytical()
