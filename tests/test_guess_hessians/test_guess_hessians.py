#!/usr/bin/env python3

from pysisyphus.calculators.XTB import XTB
from pysisyphus.helpers import geom_from_library
from pysisyphus.optimizers.guess_hessians import lindh_guess
from pysisyphus.optimizers.RFOptimizer import RFOptimizer

import matplotlib.pyplot as plt
import numpy as np


def plot_hessian(H):
    fig, ax = plt.subplots()
    im = ax.imshow(H)
    fig.colorbar(im)
    plt.show()



def test_lindh_hessian():
    geom = geom_from_library("birkholz/vitamin_c.xyz", coord_type="redund")
    geom.set_calculator(XTB())
    opt_kwargs = {
        "max_cycles": 125,
        "thresh": "gau",
        "trust_update": True,
        "hessian_init": "lindh",
        # "hessian_init": "guess",
        "trust_radius": 0.5,
        # "hessian_update": "flowchart",
        "hessian_update": "bfgs",
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()

    from pysisyphus.optimizers.RFOptimizer_org import RFOptimizer as RFOptimizer_
    geom = geom_from_library("birkholz/vitamin_c.xyz", coord_type="redund")
    geom.set_calculator(XTB())
    opt_kwargs = {
        "max_cycles": 125,
        "thresh": "gau",
        "trust_update": True,
        "trust_radius": 0.5,
    }
    opt = RFOptimizer_(geom, **opt_kwargs)
    opt.run()
    return
    # H = geom.hessian
    # np.savetxt("vc_xtb_hess", H)
    # calc = geom.calculator
    # og = calc.run_opt(geom.atoms, geom.coords)
    # import pdb; pdb.set_trace()
    H_lindh = lindh_guess(geom)
    print(np.diag(H_lindh))
    plot_hessian(H_lindh)


if __name__ == "__main__":
    test_lindh_hessian()
