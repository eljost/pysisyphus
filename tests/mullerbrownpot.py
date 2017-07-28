#!/usr/bin/env python3

import numpy as np

from AnimPlot import AnimPlot
from calculators.MullerBrownPot import MullerBrownPot
from cos.NEB import NEB
from cos.SimpleZTS import SimpleZTS
from Geometry import Geometry
from optimizers.SteepestDescent import SteepestDescent
#from optimizers.ConjugateGradient import ConjugateGradient
from optimizers.FIRE import FIRE
from optimizers.NaiveSteepestDescent import NaiveSteepestDescent

CYCLES = 50
IMAGES = 10


def get_geoms():
    min_a = np.array((-0.558, 1.442, 0)) # Minimum A
    min_b = np.array((0.6215, 0.02838, 0)) # Minimum B
    min_c = np.array((-0.05, 0.467, 0)) # Minimum C
    saddle_a = np.array((-0.822, 0.624, 0)) # Saddle point A
    coords = (min_b, min_c, saddle_a, min_a)
    #coords = (min_b, min_c)
    """
    t1 = np.array((0.4879, 0.4509, 0))
    t2 = np.array((-0.2621, 0.1027, 0))
    coords = (t1, t2)
    """
    atoms = ("H", )
    geoms = [Geometry(atoms, c) for c in coords]
    return geoms


def run_cos_opt(cos):
    cos.interpolate(IMAGES)
    for img in cos.images:
        img.set_calculator(MullerBrownPot())

    #sd = NaiveSteepestDescent(cos,
    #sd = SteepestDescent(cos,
    opt = FIRE(cos,
                         max_cycles=CYCLES,
                         #max_step=0.5,
                         alpha=0.05)
    opt.run()
    xlim = (-1.75, 1.25)
    ylim = (-0.5, 2.25)
    levels=(-150, -15, 40)
    ap = AnimPlot(MullerBrownPot(), opt, xlim=xlim, ylim=ylim, levels=levels)
    ap.animate()


if __name__ == "__main__":
    geoms = get_geoms()
    neb = NEB(geoms)
    szts = SimpleZTS(geoms, param="equal")
    szts = SimpleZTS(geoms, param="energy")

    #run_cos_opt(neb)
    print()
    run_cos_opt(szts)
