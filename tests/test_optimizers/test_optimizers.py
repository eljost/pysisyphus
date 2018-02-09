#!/usr/bin/env python3

import matplotlib.pyplot as plt

from pysisyphus.calculators.AnaPot4 import AnaPot4
from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


def test_rfoptimizer():
    atoms = ("X", )
    # http://www.applied-mathematics.net/optimization/optimizationIntro.html
    coords = (0.7, -3.3, 0)
    geom = Geometry(atoms, coords)
    ap4 = AnaPot4()
    geom.set_calculator(ap4)
    opt = RFOptimizer(geom)
    ap4.plot()
    plt.show()


if __name__ == "__main__":
    test_rfoptimizer()
