#!/usr/bin/env python3

import copy

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pysisyphus.plotters.AnimPlot import AnimPlot
from pysisyphus.calculators.AnaPot3 import AnaPot3
from pysisyphus.cos.NEB import NEB
from pysisyphus.cos.SimpleZTS import SimpleZTS
from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.BFGS import BFGS
from pysisyphus.optimizers.FIRE import FIRE
from pysisyphus.optimizers.SteepestDescent import SteepestDescent
from pysisyphus.optimizers.NaiveSteepestDescent import NaiveSteepestDescent

KWARGS = {
    "images": 10,
    "max_cycles": 50,
    "convergence": {
        "max_force_thresh": 1.9e-2,
        "rms_force_thresh": 9.e-3,
        "max_step_thresh": 1.0e-2,
        "rms_step_thresh": 4e-3,
    },
    "dump": False,
}


def get_geoms():
    initial = np.array((-0.5, 0.5, 0))
    final = np.array((0.5, 0.5, 0))
    coords = (initial, final)
    atoms = ("H")
    geoms = [Geometry(atoms, c) for c in coords]
    return geoms


def run_cos_opt(cos, Opt, images, **kwargs):
    cos.interpolate(images)
    opt = Opt(cos, **kwargs)
    for img in cos.images:
        img.set_calculator(AnaPot3())
    opt.run()

    return opt


def animate(opt):
    xlim = (-1.5, 1.5)
    ylim = (-0.5, 1.5)
    levels = (-.5, 2, 30)
    interval = 250
    ap = AnimPlot(
            AnaPot3(), opt, xlim=xlim, ylim=ylim,
            levels=levels, interval=interval
    )
    ap.animate()


@pytest.mark.sd
def test_steepest_descent_neb():
    kwargs = copy.copy(KWARGS)
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 24)

    return opt


@pytest.mark.sd
def test_steepest_descent_neb_more_images():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 20
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 29)

    return opt


@pytest.mark.fire
def test_fire_neb():
    kwargs = copy.copy(KWARGS)
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, FIRE, **kwargs)
    
    assert(opt.is_converged)
    assert(opt.cur_cycle == 18)

    return opt


@pytest.mark.bfgs
@pytest.mark.skip("Doesn't work at all!")
def test_bfgs_neb():
    kwargs = copy.copy(KWARGS)
    #kwargs["max_cycles"] = 18
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, BFGS, **kwargs)

    assert(opt.is_converged)

    return opt


@pytest.mark.bfgs
@pytest.mark.skip("Doesn't work at all!")
def test_bfgs_neb_more_images():
    kwargs = copy.copy(KWARGS)
    kwargs["max_cycles"] = 18
    kwargs["images"] = 10
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, BFGS, **kwargs)

    assert(opt.is_converged)

    return opt


def test_equal_szts():
    kwargs = copy.copy(KWARGS)
    convergence = {
        "max_force_thresh": 0.06,
    }
    kwargs["convergence"] = convergence
    szts_equal = SimpleZTS(get_geoms(), param="equal")
    opt = run_cos_opt(szts_equal, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 20)

    return opt


def test_equal_szts_more_images():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 20
    convergence = {
        "max_force_thresh": 0.05,
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
        "max_force_thresh": 9.2e-2,
    }
    kwargs["convergence"] = convergence
    szts_energy = SimpleZTS(get_geoms(), param="energy")
    opt = run_cos_opt(szts_energy, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 20)

    return opt


def test_energy_szts_more_images():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 20
    convergence = {
        "max_force_thresh": 0.1,
    }
    kwargs["convergence"] = convergence
    szts_energy = SimpleZTS(get_geoms(), param="energy")
    opt = run_cos_opt(szts_energy, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 20)

    return opt


if __name__ == "__main__":
    # Steepest Descent
    opt = test_steepest_descent_neb()
    #opt = test_steepest_descent_neb_more_images()

    # FIRE
    #opt = test_fire_neb()

    # BFGS fails completely here!
    #opt = test_bfgs_neb()
    #opt = test_bfgs_neb_more_images()

    # SimpleZTS
    #opt = test_equal_szts()
    #opt = test_equal_szts_more_images()
    #opt = test_energy_szts()
    #opt = test_energy_szts_more_images()

    ap = animate(opt)
    plt.show()
