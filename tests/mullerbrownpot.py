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
IMAGES = 5


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

    #opt = NaiveSteepestDescent(cos,
    #opt = SteepestDescent(cos,
    opt = FIRE(cos,
                         max_cycles=CYCLES)
                         #max_step=0.005,
                         #dt_max=0.2)
    opt.run()
    xlim = (-1.75, 1.25)
    ylim = (-0.5, 2.25)
    levels=(-150, -15, 40)
    ap = AnimPlot(MullerBrownPot(), opt, xlim=xlim, ylim=ylim, levels=levels)
    ap.animate()


if __name__ == "__main__":
    """NEB doesn't converge with the default thresholds at all but the
    SimpleZTS methods converge fine."""
    neb = NEB(get_geoms())
    run_cos_opt(neb)
    print()

    szts_equal = SimpleZTS(get_geoms(), param="equal")
    run_cos_opt(szts_equal)
    print()

    szts_energy = SimpleZTS(get_geoms(), param="energy")
    run_cos_opt(szts_energy)
