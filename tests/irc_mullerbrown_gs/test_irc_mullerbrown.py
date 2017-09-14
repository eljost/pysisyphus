#!/usr/bin/env python3

import os
import pathlib

import numpy as np

from pysisyphus.calculators.MullerBrownSympyPot2D import MullerBrownSympyPot2D
from pysisyphus.Geometry import Geometry
from pysisyphus.irc.GonzalesSchlegel import GonzalesSchlegel
from pysisyphus.irc.DampedVelocityVerlet import DampedVelocityVerlet

def test_mullerbrown_gs_irc():
    this_dir = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    atoms = ("H", )
    # Original ts coordinates
    calc, ts_coords = (MullerBrownSympyPot2D(), np.array((-0.822, 0.624)))
    geometry = Geometry(atoms, ts_coords)
    geometry.set_calculator(calc)

    gs_irc = GonzalesSchlegel(geometry, step_length=0.3, max_steps=6)
    gs_irc.run()
    assert(gs_irc.forward_step == 5)
    assert(gs_irc.backward_step == 4)

    return gs_irc

def test_mullerbrown_dvv_irc():
    this_dir = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
    atoms = ("H", )
    # Original ts coordinates
    calc, ts_coords = (MullerBrownSympyPot2D(), np.array((-0.822, 0.624)))
    geometry = Geometry(atoms, ts_coords)
    geometry.set_calculator(calc)

    #dvv = DampedVelocityVerlet(geometry, v0=1, energy_lowering=1e-2, backward=True)
    dvv = DampedVelocityVerlet(geometry, v0=0.04, max_steps=10)
    dvv.run()

    return dvv

if __name__ == "__main__":
    #gs_irc = test_mullerbrown_gs_irc()
    #gs_irc.show2d()
    dvv = test_mullerbrown_dvv_irc()
    dvv.show2d()
