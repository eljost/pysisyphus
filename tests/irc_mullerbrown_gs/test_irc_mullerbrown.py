#!/usr/bin/env python3

import os
import pathlib

import numpy as np

from pysisyphus.calculators.MullerBrownSympyPot2D import MullerBrownSympyPot2D
from pysisyphus.Geometry import Geometry
from pysisyphus.irc.GonzalesSchlegel import GonzalesSchlegel

def test_mullerbrown_gs_irc():
    this_dir = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    atoms = ("H", )
    # Original ts coordinates
    calc, ts_coords = (MullerBrownSympyPot2D(), np.array((-0.822, 0.624)))
    geometry = Geometry(atoms, ts_coords)
    geometry.set_calculator(calc)

    gs_irc = GonzalesSchlegel(geometry, max_step=0.2, energy_lowering=5e-3)
    gs_irc.run()

    return gs_irc

if __name__ == "__main__":
    gs_irc = test_mullerbrown_gs_irc()
    gs_irc.show2d()
