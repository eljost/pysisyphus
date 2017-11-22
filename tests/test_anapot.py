#!/usr/bin/env python3

import copy

import numpy as np
import pytest

from pysisyphus.AnimPlot import AnimPlot
from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.cos.NEB import NEB
from pysisyphus.cos.SimpleZTS import SimpleZTS
from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.BFGS import BFGS
from pysisyphus.optimizers.ConjugateGradient import ConjugateGradient
from pysisyphus.optimizers.QuickMin import QuickMin
from pysisyphus.optimizers.FIRE import FIRE
from pysisyphus.optimizers.SteepestDescent import SteepestDescent
from pysisyphus.optimizers.NaiveSteepestDescent import NaiveSteepestDescent

KWARGS = {
    "images": 5,
    "max_cycles": 100,
    "convergence": {
        "max_force_thresh": 6e-3,
        "rms_force_thresh": 5e-3,
        "max_step_thresh": 6e-4,
        "rms_step_thresh": 3e-4,
    },
    "dump": False,
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


def animate(opt):
    xlim = (-2, 2.5)
    ylim = (0, 5)
    levels = (-3, 4, 80)
    ap = AnimPlot(AnaPot(), opt, xlim=xlim, ylim=ylim, levels=levels)
    ap.animate()


@pytest.mark.sd
def test_steepest_descent_neb():
    kwargs = copy.copy(KWARGS)
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 23) # k = 0.01

    return opt


@pytest.mark.sd
def test_fixed_ends_neb():
    coords = np.array(((-0.916, 1.034, 0), (1.85, 3.57, 0)))
    kwargs = copy.copy(KWARGS)
    neb = NEB(get_geoms(coords), fix_ends=True)
    convergence = {
        "max_force_thresh": 2.3e-2,
        "rms_force_thresh": 2.1e-2,
        "max_step_thresh": 8.0e-4,
        "rms_step_thresh": 5.2e-4,
    }
    kwargs["convergence"] = convergence
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 22) # k = 0.01

    return opt


@pytest.mark.sd
def test_fixed_end_climbing_neb():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 15
    convergence = {
        "max_force_thresh": 8.9e-1,
        "rms_force_thresh": 1.5e-1,
        "max_step_thresh": 5.0e-2,
        "rms_step_thresh": 8.2e-3,
    }
    kwargs["convergence"] = convergence
    neb = NEB(get_geoms(), fix_ends=True,
              k_max=5, k_min=1,
              climb=True, climb_after=5)
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 29)

    return opt


@pytest.mark.sd
def test_steepest_descent_neb_more_images():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 10
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 21) # k = 0.01

    return opt


@pytest.mark.cg
@pytest.mark.skip
def test_cg_neb():
    kwargs = copy.copy(KWARGS)
    kwargs["max_cycles"] = 25
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, ConjugateGradient, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 25) # k = 0.01

    return opt


@pytest.mark.qm
def test_qm_neb():
    kwargs = copy.copy(KWARGS)
    #kwargs["max_cycles"] = 22
    kwargs["dt"] = 0.1
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, QuickMin, **kwargs)

    #assert(opt.is_converged)
    #assert(opt.cur_cycle == 20) # k = 0.01

    return opt


@pytest.mark.fire
def test_fire_neb():
    kwargs = copy.copy(KWARGS)
    #kwargs["max_cycles"] = 23
    kwargs["dt"] = 0.05
    kwargs["dt_max"] = 0.5
    convergence = {
        "max_force_thresh": 0.0013,
        "rms_force_thresh": 0.0012,
        "max_step_thresh": 2.1e-3,
        "rms_step_thresh": 6.6e-4,
    }
    kwargs["convergence"] = convergence
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, FIRE, **kwargs)
    
    assert(opt.is_converged)
    assert(opt.cur_cycle == 30) # k = 0.01

    return opt


@pytest.mark.bfgs
def test_bfgs_neb():
    kwargs = copy.copy(KWARGS)
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, BFGS, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 23) # k = 0.01

    return opt


@pytest.mark.bfgs
def test_bfgs_neb_more_images():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 10
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, BFGS, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 32) # k = 0.01

    return opt


def test_equal_szts():
    kwargs = copy.copy(KWARGS)
    convergence = {
        "max_force_thresh": 0.1,
    }
    kwargs["convergence"] = convergence
    szts_equal = SimpleZTS(get_geoms(), param="equal")
    opt = run_cos_opt(szts_equal, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 26) # k = 0.01

    return opt


def test_equal_szts_more_images():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 10
    convergence = {
        "max_force_thresh": 0.1,
    }
    kwargs["convergence"] = convergence
    szts_equal = SimpleZTS(get_geoms(), param="equal")
    opt = run_cos_opt(szts_equal, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 28) # k = 0.01

    return opt


def test_energy_szts():
    kwargs = copy.copy(KWARGS)
    convergence = {
        "max_force_thresh": 0.1,
    }
    kwargs["convergence"] = convergence
    szts_energy = SimpleZTS(get_geoms(), param="energy")
    opt = run_cos_opt(szts_energy, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 27) # k = 0.01

    return opt


def test_energy_szts_more_images():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 10
    convergence = {
        "max_force_thresh": 0.1,
    }
    kwargs["convergence"] = convergence
    szts_energy = SimpleZTS(get_geoms(), param="energy")
    opt = run_cos_opt(szts_energy, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 28) # k = 0.01

    return opt


if __name__ == "__main__":
    # Steepest Descent
    #opt = test_steepest_descent_neb()
    #opt = test_fixed_ends_neb()
    #opt = test_fixed_end_climbing_neb()
    #opt = test_steepest_descent_neb_more_images()

    # Conjugate Gradient
    #opt = test_cg_neb()

    # QuickMin
    opt = test_qm_neb()

    # FIRE
    #opt = test_fire_neb()

    # BFGS
    #opt = test_bfgs_neb()
    #opt = test_bfgs_neb_more_images()

    # SimpleZTS
    #opt = test_equal_szts()
    #opt = test_equal_szts_more_images()
    #opt = test_energy_szts()
    #opt = test_energy_szts_more_images()

    animate(opt)
