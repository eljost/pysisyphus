#!/usr/bin/env python3

import copy

import numpy as np
from pytest import approx

from pysisyphus.AnimPlot import AnimPlot
from pysisyphus.calculators.MullerBrownPot import MullerBrownPot
from pysisyphus.cos.NEB import NEB
from pysisyphus.cos.SimpleZTS import SimpleZTS
from pysisyphus.optimizers.FIRE import FIRE
from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.SteepestDescent import SteepestDescent

KWARGS = {
    "images": 4,
    "max_cycles": 100,
    "max_step": 0.02,
    "convergence": {
        "max_force_thresh": 1.3,
        "rms_force_thresh": 0.3,
        "max_step_thresh": 5e-4,
        "rms_step_thresh": 2e-4,
    },
}


def get_geoms(coords=None):
    if coords is None:
        min_a = np.array((-0.558, 1.442, 0)) # Minimum A
        min_b = np.array((0.6215, 0.02838, 0)) # Minimum B
        min_c = np.array((-0.05, 0.467, 0)) # Minimum C
        saddle_a = np.array((-0.822, 0.624, 0)) # Saddle point A
        coords = (min_b, min_c, saddle_a, min_a)
    atoms = ("H")
    geoms = [Geometry(atoms, c) for c in coords]
    return geoms


def run_cos_opt(cos, Opt, images, **kwargs):
    cos.interpolate(images)
    opt = Opt(cos, **kwargs)
    for img in cos.images:
        img.set_calculator(MullerBrownPot())
    opt.run()

    return opt


def animate(opt):
    xlim = (-1.75, 1.25)
    ylim = (-0.5, 2.25)
    levels=(-150, -15, 40)
    ap = AnimPlot(MullerBrownPot(), opt, xlim=xlim, ylim=ylim, levels=levels)
    ap.animate()


def test_steepest_descent_neb():
    kwargs = copy.copy(KWARGS)
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 26)

    return opt


def test_climbing_neb():
    kwargs = copy.copy(KWARGS)
    neb = NEB(get_geoms(), climb=True, climb_after=15)
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    #assert(opt.is_converged)
    #assert(opt.cur_cycle == 23) # k = 0.01

    return opt

    
def test_steepest_descent_neb_more_images():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 7
    convergence = {
        "max_force_thresh": 1.8,
        "rms_force_thresh": 0.38,
        "max_step_thresh": 3e-3,
        "rms_step_thresh": 6e-4,
    }
    kwargs["convergence"] = convergence
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 19)

    return opt



def test_fire_neb():
    kwargs = copy.copy(KWARGS)
    kwargs["dt"] = 0.01
    kwargs["dt_max"] = 0.1
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, FIRE, **kwargs)
 
    assert(opt.is_converged)
    assert(opt.cur_cycle == 17)

    return opt


def test_equal_szts():
    kwargs = copy.copy(KWARGS)
    convergence = {
        "rms_force_thresh": 2.4,
    }
    kwargs["convergence"] = convergence
    szts_equal = SimpleZTS(get_geoms(), param="equal")
    opt = run_cos_opt(szts_equal, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 18)

    return opt


def test_equal_szts_straight():
    coords = np.array(((0, -0.25, 0), (0, 1.5, 0)))
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 10
    kwargs["max_step"] = 0.04
    convergence = {
        "rms_force_thresh": 2.4,
    }
    kwargs["convergence"] = convergence
    szts_equal = SimpleZTS(get_geoms(coords), param="equal")
    opt = run_cos_opt(szts_equal, SteepestDescent, **kwargs)

    return opt


def test_equal_szts_more_images():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 7
    convergence = {
        "rms_force_thresh": 2.4,
    }
    kwargs["convergence"] = convergence
    szts_equal = SimpleZTS(get_geoms(), param="equal")
    opt = run_cos_opt(szts_equal, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 22)

    return opt


def test_energy_szts():
    kwargs = copy.copy(KWARGS)
    convergence = {
        "rms_force_thresh": 2.8,
    }
    kwargs["convergence"] = convergence
    szts_energy = SimpleZTS(get_geoms(), param="energy")
    opt = run_cos_opt(szts_energy, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 16)

    return opt


def test_energy_szts_more_images():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 10
    convergence = {
        "rms_force_thresh": 1.7,
    }
    kwargs["convergence"] = convergence
    szts_energy = SimpleZTS(get_geoms(), param="energy")
    opt = run_cos_opt(szts_energy, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 23)

    return opt

if __name__ == "__main__":
    # Steepest Descent
    #opt = test_steepest_descent_neb()
    opt = test_climbing_neb()
    #opt = test_steepest_descent_neb_more_images()

    # FIRE
    #opt = test_fire_neb()

    # SimpleZTS
    #opt = test_equal_szts()
    #opt = test_equal_szts_straight()
    #opt = test_equal_szts_more_images()
    #opt = test_energy_szts()
    #opt = test_energy_szts_more_images()

    animate(opt)
