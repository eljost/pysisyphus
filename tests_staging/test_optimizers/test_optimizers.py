#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.AnaPot4 import AnaPot4
from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.SteepestDescent import SteepestDescent
from pysisyphus.optimizers.BFGS import BFGS
from pysisyphus.plotters.RFOPlotter import RFOPlotter


def test_rfoptimizer():
    kwargs = {
        #"max_cycles": 10,
        "trust_radius": 0.5,
        "max_step": 0.5,
    }
    atoms = ("X", )
    # http://www.applied-mathematics.net/optimization/optimizationIntro.html
    coords = (0.7, -3.3, 0)
    geom = Geometry(atoms, coords)
    ap4 = AnaPot4()
    geom.set_calculator(ap4)
    opt = RFOptimizer(geom, **kwargs)
    #opt = SteepestDescent(geom, **kwargs)
    #opt = BFGS(geom, max_step=1.0)
    opt.run()
    rfop = RFOPlotter(ap4, opt, save=True, title=False)
    rfop.plot()
    plt.show()


if __name__ == "__main__":
    test_rfoptimizer()
