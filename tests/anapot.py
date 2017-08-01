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
IMAGES = 5

def get_geoms():
    initial = np.array((-1.05274, 1.02776, 0))
    final = np.array((1.94101, 3.85427, 0))
    #initial = np.array((0.3, 2.46, 0))
    #final = np.array((0.8833, 0.5238, 0))
    #initial = np.array((0.3167, 2.2381, 0.0))
    #final = np.array((0.8583, 0.7381, 0.0))
    atoms = ("H")
    geoms = [Geometry(atoms, coords) for coords in (initial, final)]
    return geoms


def run_cos_opt(cos, **kwargs):
    cos.interpolate(IMAGES)
    for img in cos.images:
        img.set_calculator(AnaPot())

    opt = FIRE(cos, **kwargs)
    #opt = SteepestDescent(cos, **kwargs)
    opt.run()

    xlim = (-2, 2.5)
    ylim = (0, 5)
    levels = (-6, 6, 50)
    ap = AnimPlot(AnaPot(), opt, xlim=xlim, ylim=ylim, levels=levels)
    ap.animate()


if __name__ == "__main__":
    kwargs = {
        "max_cycles": CYCLES,
        #"max_step": 0.5,
        #"dt_max": 0.3,
    }

    neb = NEB(get_geoms())
    run_cos_opt(neb, **kwargs)
    print()

    #szts_equal = SimpleZTS(get_geoms(), param="equal")
    #run_cos_opt(szts_equal, **kwargs)
    #print()

    #szts_energy = SimpleZTS(get_geoms(), param="energy")
    #run_cos_opt(szts_energy)
