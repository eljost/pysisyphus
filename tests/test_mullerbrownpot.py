#!/usr/bin/env python3

import numpy as np

from pysisyphus.AnimPlot import AnimPlot
from calculators.MullerBrownPot import MullerBrownPot
from pysisyphus.cos.NEB import NEB
from pysisyphus.cos.SimpleZTS import SimpleZTS
from pysisyphus.optimizers.FIRE import FIRE
from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.SteepestDescent import SteepestDescent

CYCLES = 50
IMAGES = 5
KWARGS = {
    "max_cycles": CYCLES,
}


def get_geoms():
    min_a = np.array((-0.558, 1.442, 0)) # Minimum A
    min_b = np.array((0.6215, 0.02838, 0)) # Minimum B
    min_c = np.array((-0.05, 0.467, 0)) # Minimum C
    saddle_a = np.array((-0.822, 0.624, 0)) # Saddle point A
    coords = (min_b, min_c, saddle_a, min_a)
    atoms = ("H")
    geoms = [Geometry(atoms, c) for c in coords]
    return geoms


def run_cos_opt(cos, Opt, **kwargs):
    cos.interpolate(IMAGES)
    opt = Opt(cos, **kwargs)
    for img in cos.images:
        img.set_calculator(MullerBrownPot())
    opt.run()

    return opt


def animate(opt):
    xlim = (-1.75, 1.25)
    ylim = (-0.5, 2.25)
    levels=(-150, -15, 40)
    ap = AnimPlot(AnaPot(), opt, xlim=xlim, ylim=ylim, levels=levels)
    ap.animate()


def test_steepest_descent_neb():
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, SteepestDescent, **KWARGS)
    return opt


def test_fire_neb():
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, FIRE, **KWARGS)
    return opt


def test_equal_szts():
    szts_equal = SimpleZTS(get_geoms(), param="equal")
    opt = run_cos_opt(szts_equal, SteepestDescent, **KWARGS)
    return opt


def test_energy_szts():
    szts_energy = SimpleZTS(get_geoms(), param="energy")
    opt = run_cos_opt(szts_energy, SteepestDescent, **KWARGS)
    return opt

if __name__ == "__main__":
    opt = test_energy_szts()
    animate(opt)
