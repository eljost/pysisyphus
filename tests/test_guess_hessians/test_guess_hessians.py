#!/usr/bin/env python3

from pysisyphus.calculators.XTB import XTB
from pysisyphus.helpers import (geom_from_xyz_file,
                                geom_from_library,
                                do_final_hessian,
                               )
from pysisyphus.optimizers.guess_hessians import (lindh_guess, fischer_guess,
                                                  simple_guess, swart_guess,
                                                  ts_hessian)
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer

import matplotlib.pyplot as plt
import numpy as np


def plot_hessian(H, title=""):
    fig, ax = plt.subplots()
    im = ax.imshow(H, vmin=0, vmax=1)
    fig.colorbar(im)
    fig.suptitle(title)


def test_guess_hessians():
    geom = geom_from_library("birkholz/vitamin_c.xyz", coord_type="redund")

    H_lindh = lindh_guess(geom)
    print("Lindh")
    print(np.diag(H_lindh))
    plot_hessian(H_lindh, "Lindh")

    H_fischer = fischer_guess(geom)
    print("Fischer")
    print(np.diag(H_fischer))
    plot_hessian(H_fischer, "Fischer")

    H_simple = simple_guess(geom)
    print("Simple")
    print(np.diag(H_simple))
    plot_hessian(H_simple, "Simple 0.5/0.2/0.1")

    H_swart = swart_guess(geom)
    print("Swart")
    print(np.diag(H_swart))
    plot_hessian(H_swart, "Swart")

    plt.show()


def test_ts_hessian():
    H = np.diag((1, 0.5, 0.25))
    tsh = ts_hessian(H, (0, 1))
    w, v = np.linalg.eigh(tsh)
    np.testing.assert_allclose(w, (-0.445194, 0.070194, 0.25), atol=1e-7)


def test_ts_hessian_opt():
    geom = geom_from_xyz_file("/scratch/Code/parsezmat/01_hcn.xyz", coord_type="redund",)
    # geom = geom_from_library("hcn_bent_stretched.xyz", coord_type="redund")
    # geom = geom_from_library("hcn_bent_stretched.xyz")
    geom.set_calculator(XTB(pal=4))

    opt_kwargs = {
        # "hessian_init": "fischer",
        "dump": True,
        "trust_radius": 0.3,
        "trust_max": 0.3,
        "rx_coords": ((2, 1, 0), )
    }
    opt = RSPRFOptimizer(geom, **opt_kwargs)
    opt.run()
    do_final_hessian(geom)


if __name__ == "__main__":
    # test_guess_hessians()
    # test_ts_hessian()
    test_ts_hessian_opt()
