#!/usr/bin/env python3

import numpy as np

from AnimPlot import AnimPlot
from calculators.AnaPot import AnaPot
from cos.NEB import NEB
from cos.SimpleZTS import SimpleZTS
from Geometry import Geometry
from optimizers.SteepestDescent import SteepestDescent

CYCLES = 25
IMAGES = 15

def get_geoms():
    initial = np.array((-1.05274, 1.02776, 0))
    final = np.array((1.94101, 3.85427, 0))
    atoms = ("H", "H")
    geoms = [Geometry(atoms, coords) for coords in (initial, final)]
    return geoms


def run_cos_opt(cos_class, reparametrize=False, animate=False):
    geoms = get_geoms()
    cos = cos_class(geoms)
    cos.interpolate(IMAGES)
    for img in cos.images:
        img.set_calculator(AnaPot())

    #sd = NaiveSteepestDescent(cos,
    sd = SteepestDescent(cos,
                         max_cycles=CYCLES,
                         max_force_thresh=0.05,
                         rms_force_thresh=0.01,
                         alpha=-0.05)
    if reparametrize:
        sd.run(reparam=cos.reparametrize)
    else:
        sd.run()

    if animate:
        xlim = (-2, 2.5)
        ylim = (0, 5)
        levels = (-4, 8, 20)
        ap = AnimPlot(AnaPot(), sd, xlim=xlim, ylim=ylim, levels=levels)
        ap.animate()


if __name__ == "__main__":
    #run_cos_opt(NEB)
    print()
    run_cos_opt(SimpleZTS, reparametrize=True, animate=True)

