#!/usr/bin/env python3

import numpy as np

from pysisyphus.calculators.MullerBrownSympyPot2D import MullerBrownSympyPot2D
from pysisyphus.Geometry import Geometry
from pysisyphus.irc.GonzalesSchlegel import GonzalesSchlegel

def run():
    atoms = ("H", )
    # Original ts coordinates
    calc, ts_coords = (MullerBrownSympyPot2D(), np.array((-0.822, 0.624)))
    geometry = Geometry(atoms, ts_coords)
    geometry.set_calculator(calc)

    GS = GonzalesSchlegel(geometry, max_step=0.2)
    GS.run()
    GS.show2d()

if __name__ == "__main__":
    run()
