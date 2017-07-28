#!/usr/bin/env python3

import numpy as np

from AnimPlot import AnimPlot
from calculators.AnaPot import AnaPot
from cos.NEB import NEB
from cos.SimpleZTS import SimpleZTS
from optimizers.FIRE import FIRE
from Geometry import Geometry
from optimizers.SteepestDescent import SteepestDescent

CYCLES = 50
IMAGES = 15

def get_geoms():
    initial = np.array((-1.05274, 1.02776, 0))
    final = np.array((1.94101, 3.85427, 0))
    atoms = ("H", "H")
    geoms = [Geometry(atoms, coords) for coords in (initial, final)]
    return geoms


def run_cos_opt(cos):
    cos.interpolate(IMAGES)
    for img in cos.images:
        img.set_calculator(AnaPot())

    #opt = SteepestDescent(cos,
    opt = FIRE(cos,
                         max_cycles=CYCLES,
                         max_step = 0.05,
                         #max_force_thresh=0.05,
                         #rms_force_thresh=0.01,
                         alpha=0.1)
    opt.run()

    xlim = (-2, 2.5)
    ylim = (0, 5)
    levels = (-4, 8, 20)
    ap = AnimPlot(AnaPot(), opt, xlim=xlim, ylim=ylim, levels=levels)
    ap.animate()


if __name__ == "__main__":
    geoms = get_geoms()
    neb = NEB(geoms)
    szts = SimpleZTS(geoms, param="equal")
    szts = SimpleZTS(geoms, param="energy")

    run_cos_opt(neb)
    print()
    run_cos_opt(szts)
