#!/usr/bin/env python3

import copy

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pysisyphus.plotters.AnimPlot import AnimPlot
from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.cos.NEB import NEB
from pysisyphus.cos.SimpleZTS import SimpleZTS
from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.BFGS import BFGS
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.ConjugateGradient import ConjugateGradient
from pysisyphus.optimizers.QuickMin import QuickMin
from pysisyphus.optimizers.FIRE import FIRE
from pysisyphus.optimizers.SteepestDescent import SteepestDescent
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.SciPyOptimizer import SciPyOptimizer
from pysisyphus.optimizers.closures import modified_broyden_closure


KWARGS = {
    "images": 5,
    "max_cycles": 50,
    "rms_force": 5e-3,
    "dump": False,
    "align": False,
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
    return ap


def animate_bare(opt):
    xlim = (-2, 2.5)
    ylim = (0, 5)
    levels = (-3, 4, 80)
    ap = AnimPlot(AnaPot(), opt, xlim=xlim, ylim=ylim, levels=levels,
                  energy_profile=False, colorbar=False, figsize=(8, 6),
                  save=False, title=False,
    )
    ap.animate()
    return ap

@pytest.mark.sd
def test_steepest_descent_neb():
    kwargs = copy.copy(KWARGS)
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 29)  # k = 0.01

    return opt


@pytest.mark.sd
def test_spline_hei():
    kwargs = copy.copy(KWARGS)
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 29)  # k = 0.01

    hei_coords, hei_energy, hei_tangent = neb.get_splined_hei()
    # ap = animate(opt)
    # hei = hei_coords[:2]
    # heit = hei_tangent[:2]
    # ap.ax.scatter(*hei, s=15, c="k")
    # ap.ax.quiver(*hei, *heit)
    # plt.show()
    print(f"Interpolated HEI: {hei_coords}")
    print(f"HEI tangent: {hei_tangent}")

    from pysisyphus.tsoptimizers.dimer import dimer_method
    calc_getter = AnaPot
    ts_geom = neb.images[0].copy()
    ts_geom.coords = hei_coords
    ts_geom.set_calculator(calc_getter())
    geoms = [ts_geom, ]
    dimer_kwargs = {
        "ana_2dpot": True,
        "restrict_step": "max",
        "N_init": hei_tangent,
        "calc_getter": calc_getter,
    }
    dimer_cycles = dimer_method(geoms, **dimer_kwargs)
    last_cycle = dimer_cycles[-1]
    ts_coords = last_cycle.trans_coords[1]
    print(f"Optimized TS coord: {ts_coords}")

    return opt


"""
@pytest.mark.sd
def test_steepest_descent_var_springs_neb():
    kwargs = copy.copy(KWARGS)
    kwargs["climb"] = True
    kwargs["climb_rms"] = 0.01
    convergence = {
        "max_force_thresh": 8e-3,
        "rms_force_thresh": 3.3e-3,
        "max_step_thresh": 3.2e-3,
        "rms_step_thresh": 1.6e-2,
    }
    kwargs["convergence"] = convergence
    var_springs = True
    neb = NEB(get_geoms(), variable_springs=var_springs, fix_ends=True)
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 22) # k = 0.01

    return opt
"""


@pytest.mark.sd
def test_steepest_descent_neb_more_images():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 10
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 31)  # k = 0.01

    return opt


@pytest.mark.sd
def test_fix_first_neb():
    # First image is fixed at a non equilibrium geometry.
    coords = np.array(((-0.916, 1.034, 0), (1.94101, 3.85427, 0)))
    kwargs = copy.copy(KWARGS)
    neb = NEB(get_geoms(coords), fix_first=True)
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 29)  # k = 0.01

    return opt


@pytest.mark.sd
def test_fix_last_neb():
    # Last image is fixed at a non equilibrium geometry.
    coords = np.array(((-1.05274, 1.02776, 0), (1.85, 3.57, 0)))
    kwargs = copy.copy(KWARGS)
    neb = NEB(get_geoms(coords), fix_last=True)
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 28)  # k = 0.01

    return opt


@pytest.mark.sd
def test_fix_ends_neb():
    kwargs = copy.copy(KWARGS)
    neb = NEB(get_geoms(), fix_ends=True)
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 29)  # k = 0.01

    return opt


@pytest.mark.sd
def test_fix_displaced_ends_neb():
    coords = np.array(((-0.916, 1.034, 0), (1.85, 3.57, 0)))
    kwargs = copy.copy(KWARGS)
    neb = NEB(get_geoms(coords), fix_ends=True)
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 25)  # k = 0.01

    return opt


@pytest.mark.sd
def test_fix_end_immediate_climbing_neb():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 10
    kwargs["rms_force"] = 2e-2
    cos_kwargs = {
        "climb": True,
        "climb_rms": -1,
        "fix_ends": True,
    }
    neb = NEB(get_geoms(), **cos_kwargs)
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 35)

    return opt


@pytest.mark.cg
def test_cg_neb():
    kwargs = copy.copy(KWARGS)
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, ConjugateGradient, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 37)  # k = 0.01

    return opt


@pytest.mark.qm
def test_qm_neb():
    kwargs = copy.copy(KWARGS)
    kwargs["dt"] = 0.1
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, QuickMin, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 25)  # k = 0.01

    return opt


@pytest.mark.fire
def test_fire_neb():
    kwargs = copy.copy(KWARGS)
    kwargs["dt_max"] = 0.2
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, FIRE, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 38)  # k = 0.01

    return opt


@pytest.mark.fire
def test_fire_climb_neb():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 10
    kwargs["dt_max"] = 0.2
    kwargs["climb"] = True
    kwargs["climb_rms"] = 0.0065
    neb = NEB(get_geoms(), fix_ends=True)
    opt = run_cos_opt(neb, FIRE, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 39)  # k = 0.01

    return opt


@pytest.mark.bfgs
def test_bfgs_neb():
    kwargs = copy.copy(KWARGS)
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, BFGS, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 35)  # k = 0.01

    return opt


@pytest.mark.bfgs
def test_bfgs_neb_more_images():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 10
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, BFGS, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 36)  # k = 0.01

    return opt


@pytest.mark.skip
def test_fix_end_climbing_bfgs_neb():
    kwargs = copy.copy(KWARGS)
    # kwargs["bt_disable"] = True
    neb = NEB(get_geoms(), climb=True, fix_ends=True)
    opt = run_cos_opt(neb, BFGS, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 32)

    return opt


@pytest.mark.skip
def test_scipy_bfgs_neb():
    """It seems like the maximum step size can't be set
    as of SciPy 1.0.0."""
    kwargs = copy.copy(KWARGS)
    kwargs["method"] = "BFGS"
    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, SciPyOptimizer, **kwargs)

    return opt


@pytest.mark.lbfgs
def test_lbfgs_neb():
    kwargs = copy.copy(KWARGS)
    # kwargs["bt_disable"] = True
    # kwargs["rms_force"] = 1e-5
    # kwargs["alpha"] = 0.5
    kwargs["max_cycles"] = 25

    neb = NEB(get_geoms())
    opt = run_cos_opt(neb, LBFGS, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 22)  # k = 0.01

    return opt


@pytest.mark.lbfgs
def test_lbfgs_mod_neb():
    from pysisyphus.optimizers.LBFGS_mod import LBFGS as LBFGSm
    kwargs = copy.copy(KWARGS)
    # kwargs["bt_disable"] = True
    kwargs["rms_force"] = 1e-5
    # kwargs["alpha"] = 0.5
    kwargs["max_cycles"] = 25
    kwargs["images"] = 3
    # kwargs["bt_disable"] = True

    coords = (
        (-1.05274, 1.02776, 0),
        # (0.2166, 1.2619, 0),
        (0.625, 1.476, 0),
        # (1.4, 2.45, 0),
        (1.94101, 3.85427, 0),
    )

    neb = NEB(get_geoms(coords), fix_ends=True)
    # opt = run_cos_opt(neb, LBFGSm, **kwargs)
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    # assert(opt.is_converged)
    # assert(opt.cur_cycle == 22)  # k = 0.01

    return opt


@pytest.mark.modified_broyden
def test_modified_broyden():
    pot = AnaPot()
    geom = pot.get_geom((-0.8333, 2, 0))
    force_getter = lambda x: geom.get_energy_and_forces_at(x)["forces"]
    def restrict_step(step, max_len=0.5):
        norm = np.linalg.norm(step)
        if norm > max_len:
            step = max_len * step / norm
        return step
    mod_broyden = modified_broyden_closure(force_getter, restrict_step=restrict_step)
    coords = [geom.coords.copy(), ]
    for i in range(50):
        step, forces = mod_broyden(geom.coords)
        geom.coords += step
        coords.append(geom.coords.copy())
        norm = np.linalg.norm(forces)
        print(f"Cycle {i:02d}: norm={norm:.4e}")
        if norm < 1e-12:
            print("Converged")
            break
    pot.plot()
    ax = pot.ax
    coords = np.array(coords)
    xs, ys = coords.T[:2]
    ax.plot(xs, ys, "ro-")
    plt.show()


def test_rfo_optimizer():
    pot = AnaPot()
    # geom = pot.get_geom((-0.8333, 2, 0))
    geom = pot.get_geom((-1, 3, 0))
    # geom = pot.get_geom((0, 3, 0))
    opt = RFOptimizer(geom, thresh="gau_tight")
    opt.run()
    coords = np.array(opt.coords)
    pot.plot()
    ax = pot.ax
    ax.plot(coords[:,0], coords[:,1], "ro-")
    plt.show()


@pytest.mark.skip
def test_equal_szts():
    kwargs = copy.copy(KWARGS)
    szts_equal = SimpleZTS(get_geoms(), param="equal")
    opt = run_cos_opt(szts_equal, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 44)  # k = 0.01

    return opt


@pytest.mark.skip
def test_equal_szts_more_images():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 10
    szts_equal = SimpleZTS(get_geoms(), param="equal")
    opt = run_cos_opt(szts_equal, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 48)  # k = 0.01

    return opt


@pytest.mark.skip
def test_energy_szts():
    kwargs = copy.copy(KWARGS)
    szts_energy = SimpleZTS(get_geoms(), param="energy")
    opt = run_cos_opt(szts_energy, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 38)  # k = 0.01

    return opt


@pytest.mark.skip
def test_energy_szts_more_images():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 10
    szts_energy = SimpleZTS(get_geoms(), param="energy")
    opt = run_cos_opt(szts_energy, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 49)  # k = 0.01

    return opt


if __name__ == "__main__":
    # Steepest Descent
    # opt = test_steepest_descent_neb()
    # opt = test_steepest_descent_var_springs_neb()
    # opt = test_steepest_descent_neb_more_images()
    # opt = test_fix_first_neb()
    # opt = test_fix_last_neb()
    # opt = test_fix_ends_neb()
    # opt = test_fix_displaced_ends_neb()
    # Steepest descent + climbing Image
    # opt = test_fix_end_immediate_climbing_neb()

    # Conjugate Gradient
    # opt = test_cg_neb()

    # QuickMin
    # opt = test_qm_neb()

    # FIRE
    # opt = test_fire_neb()
    # opt = test_fire_climb_neb()

    # BFGS
    # opt = test_bfgs_neb()
    # opt = test_bfgs_neb_more_images()
    #opt = test_scipy_bfgs_neb()
    # BFGS + climbing Image
    # opt = test_fix_end_climbing_bfgs_neb()

    # LBFGS
    # opt = test_lbfgs_neb()
    # opt = test_lbfgs_mod_neb()

    # Modified broyden
    # test_modified_broyden()

    # Rational function optimization
    test_rfo_optimizer()

    # SimpleZTS
    # opt = test_equal_szts()
    # opt = test_equal_szts_more_images()
    # opt = test_energy_szts()
    # opt = test_energy_szts_more_images()

    # opt = test_scipy_bfgs_neb()

    # NEB with Dimer
    # test_spline_hei()

    # ap = animate(opt)
    # ap = animate_bare(opt)
    # ap.as_html5("anim.html")
    # plt.show()
