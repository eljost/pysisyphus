#!/usr/bin/env python3

import copy

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pysisyphus.plotters.AnimPlot import AnimPlot
from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.MullerBrownPot import MullerBrownPot
from pysisyphus.cos.AdaptiveNEB import AdaptiveNEB
from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.SteepestDescent import SteepestDescent
from pysisyphus.optimizers.ConjugateGradient import ConjugateGradient

KWARGS = {
    "images": 3,
    #"max_cycles": 25,
    "dump": False,
    "adapt_thresh": 0.1,
}


def get_geoms(coords=None):
    if coords is None:
        initial = np.array((-1.05274, 1.02776, 0))
        final = np.array((1.94101, 3.85427, 0))
        coords = (initial, final)
    atoms = ("H")
    geoms = [Geometry(atoms, c) for c in coords]
    return geoms


def run_cos_opt(cos, Opt, images, **kwargs):
    cos.interpolate(images)
    opt = Opt(cos, **kwargs)
    for img in cos.images:
        img.set_calculator(AnaPot())
    opt.run()

    return opt


def animate_bare(opt):
    xlim = (-2, 2.5)
    ylim = (0, 5)
    levels = (-3, 4, 80)
    ap = AnimPlot(AnaPot(), opt, xlim=xlim, ylim=ylim, levels=levels,
                  energy_profile=False, colorbar=False, figsize=(8, 6))
    ap.animate()
    return ap


def test_steepest_descent_aneb():
    kwargs = copy.copy(KWARGS)
    # keep_hei is a bad idea. needs 44 cycles
    # aneb = AdaptiveNEB(get_geoms(), keep_hei=False)
    aneb = AdaptiveNEB(get_geoms())
    opt = run_cos_opt(aneb, SteepestDescent, **kwargs)

    assert opt.is_converged
    assert opt.cur_cycle == 35


    return opt


def get_muller_brown_geoms(keys=("B","C","TSA","A")):
    coords_dict = {
        "A": (-0.558, 1.442, 0), # Minimum A
        "B": (0.6215, 0.02838, 0), # Minimum B
        "C": (-0.05, 0.467, 0), # Minimum C
        "AC": (-0.57, 0.8, 0), # Between A and C
        "TSA": (-0.822, 0.624, 0) # Saddle point A
    }
    coords = [np.array(coords_dict[k]) for k in keys]
    atoms = ("H")
    geoms = [Geometry(atoms, c) for c in coords]
    return geoms


def animate_mueller_brown(opt):
    xlim = (-1.75, 1.25)
    ylim = (-0.5, 2.25)
    levels=(-150, -15, 40)
    ap = AnimPlot(MullerBrownPot(), opt, xlim=xlim, ylim=ylim, levels=levels)
    ap.animate()


def test_mueller_brown_steepest_descent_aneb():
    kwargs = copy.copy(KWARGS)
    # Needs 49 cycles with keep_hei=False
    cos_kwargs = {
        # "climb": True,
        # "climb_rms": 20,
        "keep_hei": True,
        "adapt_between": 1,
    }
    aneb = AdaptiveNEB(get_muller_brown_geoms(("B", "A")), **cos_kwargs)

    kwargs["convergence"] = {
        "max_force_thresh": 0.03,
        "rms_force_thresh": 0.011,
        "max_step_thresh": 3e-5,
        "rms_step_thresh": 1e-5,
    }

    aneb.interpolate(3)
    opt = SteepestDescent(aneb, **kwargs)
    for img in aneb.images:
        img.set_calculator(MullerBrownPot())
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == 40

    animate_mueller_brown(opt)

    return opt


if __name__ == "__main__":
    # opt = test_steepest_descent_aneb()
    # ap = animate_bare(opt)
    # plt.show()
    test_mueller_brown_steepest_descent_aneb()
    plt.show()
