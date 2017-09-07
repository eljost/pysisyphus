#!/usr/bin/env python3

import numpy as np

from pysisyphus.calculators.MullerBrownSympyPot2D import MullerBrownSympyPot2D
from pysisyphus.Geometry import Geometry
from pysisyphus.irc.GonzalesSchlegel import GonzalesSchlegel

def run():
    atoms = ("H", )
    calc, ts_coords = (MullerBrownSympyPot2D(), np.array((-0.845041, 0.663752)))
    xlim = (-1.25, -.25)
    ylim = (0.5, 1.5)
    levels=(-150, -15, 40)
    geometry = Geometry(atoms, ts_coords)
    geometry.set_calculator(calc)

    GS = GonzalesSchlegel(geometry, max_step=0.35)
    GS.run()
    GS.show2d()

if __name__ == "__main__":
    run()
