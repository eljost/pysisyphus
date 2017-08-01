#!/usr/bin/env python3

import numpy as np

from pysisyphus.AnimPlot import AnimPlot
from pysisyphus.calculators.AnaPot import AnaPot
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
    initial = np.array((-1.05274, 1.02776, 0))
    final = np.array((1.94101, 3.85427, 0))
    atoms = ("H")
    geoms = [Geometry(atoms, coords) for coords in (initial, final)]
    return geoms


def run_cos_opt(cos, Opt, **kwargs):
    cos.interpolate(IMAGES)
    opt = Opt(cos, **kwargs)
    for img in cos.images:
        img.set_calculator(AnaPot())
    opt.run()

    return opt


def animate(opt):
    xlim = (-2, 2.5)
    ylim = (0, 5)
    levels = (-6, 6, 50)
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
