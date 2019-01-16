#!/usr/bin/env python3

import copy

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pysisyphus.plotters.AnimPlot import AnimPlot
from pysisyphus.calculators.MullerBrownPot import MullerBrownPot
#from pysisyphus.calculators.MullerBrownSympyPot import MullerBrownPot
from pysisyphus.cos.NEB import NEB
from pysisyphus.cos.SimpleZTS import SimpleZTS
from pysisyphus.optimizers.FIRE import FIRE
from pysisyphus.optimizers.BFGS import BFGS
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.SteepestDescent import SteepestDescent

KWARGS = {
    "images": 4,
    "max_cycles": 100,
    "max_step": 0.02,
    "convergence": {
        "max_force_thresh": 0.1,
        "rms_force_thresh": 0.02,
        "max_step_thresh": 0.005,
        "rms_step_thresh": 0.001,
    },
    "dump": False,
}


def get_geoms(keys=("B","C","TSA","A")):
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


@pytest.mark.sd
def test_steepest_descent_neb():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 4
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 56)

    return opt


@pytest.mark.sd
def test_steepest_descent_straight_neb():
    """Something is really really wrong here."""
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 10
    kwargs["max_cycles"] = 100
    convergence = {
        "max_force_thresh": 1.16,
        "rms_force_thresh": 0.27,
        "max_step_thresh": 0.021,
        "rms_step_thresh": 0.005,
    }
    kwargs["convergence"] = convergence
    neb = NEB(get_geoms(("A", "B")))
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 62)

    return opt


@pytest.mark.bfgs
def test_bfgs_straight_neb():
    """Something is really really wrong here."""
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 10
    convergence = {
        "max_force_thresh": 5.0,
        "rms_force_thresh": 1,
        "max_step_thresh": 0.002,
        "rms_step_thresh": 0.0006,
    }
    kwargs["convergence"] = convergence
    neb = NEB(get_geoms(("A", "B")))
    opt = run_cos_opt(neb, BFGS, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 45)

    return opt


@pytest.mark.lbfgs
def test_lbfgs_neb():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 3
    kwargs["fix_ends"] = True
    k_min = 1000
    k_max = k_min+10
    neb = NEB(get_geoms(("A", "B")), k_min=k_min, k_max=k_max, fix_ends=True)
    from pysisyphus.optimizers.ConjugateGradient import ConjugateGradient
    # from pysisyphus.optimizers.LBFGS_mod import LBFGS
    opt = run_cos_opt(neb, LBFGS, **kwargs)

    # assert(opt.is_converged)
    # assert(opt.cur_cycle == 45)

    return opt


@pytest.mark.sd
def test_steepest_descent_neb_more_images():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 7
    convergence = {
        "max_force_thresh": 0.6,
        "rms_force_thresh": 0.13,
        "max_step_thresh": 0.015,
        "rms_step_thresh": 0.0033,
    }
    kwargs["convergence"] = convergence
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 41)

    return opt


@pytest.mark.fire
def test_fire_neb():
    kwargs = copy.copy(KWARGS)
    kwargs["dt"] = 0.01
    kwargs["dt_max"] = 0.1
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, FIRE, **kwargs)
 
    assert(opt.is_converged)
    assert(opt.cur_cycle == 76)

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
    assert(opt.cur_cycle == 17)

    return opt


def test_equal_szts_straight():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 10
    kwargs["max_step"] = 0.04
    convergence = {
        "rms_force_thresh": 2.4,
    }
    kwargs["convergence"] = convergence
    szts_equal = SimpleZTS(get_geoms(("A", "B")), param="equal")
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
    assert(opt.cur_cycle == 21)

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
    assert(opt.cur_cycle == 15)

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
    assert(opt.cur_cycle == 22)

    return opt

if __name__ == "__main__":
    # Steepest Descent
    opt = test_steepest_descent_neb()
    #opt = test_steepest_descent_straight_neb()
    #opt = test_steepest_descent_neb_more_images()

    # opt = test_bfgs_straight_neb()

    # opt = test_lbfgs_neb()

    # FIRE
    #opt = test_fire_neb()

    # SimpleZTS
    #opt = test_equal_szts()
    #opt = test_equal_szts_straight()
    #opt = test_equal_szts_more_images()
    #opt = test_energy_szts()
    #opt = test_energy_szts_more_images()

    ap = animate(opt)
    plt.show()
