#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pysisyphus.plotters.AnimPlot import AnimPlot
from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.cos.NEB import NEB
from pysisyphus.cos.SimpleZTS import SimpleZTS
from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.ConjugateGradient import ConjugateGradient
from pysisyphus.optimizers.QuickMin import QuickMin
from pysisyphus.optimizers.FIRE import FIRE
from pysisyphus.optimizers.SteepestDescent import SteepestDescent
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.closures import modified_broyden_closure
from pysisyphus.interpolate.Interpolator import Interpolator


def get_geoms():
    initial = AnaPot.get_geom((-1.05274, 1.02776, 0))
    final = AnaPot.get_geom((1.94101, 3.85427, 0))
    geoms = (initial, final)
    return geoms


def run_cos_opt(geoms, between, calc_cls,
                cos_cls, cos_kwargs,
                opt_cls, opt_kwargs):
    interpol = Interpolator(geoms, between=between)
    images = interpol.interpolate_all()

    for image in images:
        image.set_calculator(AnaPot())

    cos = cos_cls(images, **cos_kwargs)

    opt = opt_cls(cos, **opt_kwargs)
    opt.run()

    return opt


def assert_cos_opt(opt, ref_cycle):
    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle


@pytest.mark.parametrize(
    "opt_cls, opt_kwargs_, neb_kwargs_, ref_cycle, between",
    [
        # (SteepestDescent, {}, {}, 30, 5),
        # (SteepestDescent, {}, {}, 32, 10),
        # (ConjugateGradient, {}, {}, 40, 5),
        # (QuickMin, {"dt": 0.1,}, {}, 27, 5),
        # (FIRE, {"dt_max": 0.2,}, {}, 42, 5),
        # (LBFGS, {}, {}, 27, 5),
])
def test_anapot_neb(opt_cls, opt_kwargs_, neb_kwargs_, ref_cycle, between):
    geoms = get_geoms()

    neb_kwargs = {
        "fix_ends": True,
        "k_min": 0.01,
    }
    neb_kwargs.update(neb_kwargs_)

    opt_kwargs = {
    }
    opt_kwargs.update(opt_kwargs_)

    opt = run_cos_opt(geoms, between, AnaPot,
                      NEB, neb_kwargs,
                      opt_cls, opt_kwargs)
    assert_cos_opt(opt, ref_cycle)

    # ap = animate(opt)
    # plt.show()



@pytest.mark.parametrize(
    "between, param, ref_cycle",
    [
        (5, "equal", 49),
        (10, "equal", 49),
        (5, "energy", 41),
        (10, "energy", 46),
])
def test_anapot_szts(between, param, ref_cycle):
    geoms = get_geoms()

    szts_kwargs = {
        "fix_ends": True,
        "param": param,
    }

    opt_cls = SteepestDescent
    opt_kwargs = {
        "max_cycles": 100,
    }
    opt = run_cos_opt(geoms, between, AnaPot,
                      SimpleZTS, szts_kwargs,
                      SteepestDescent, opt_kwargs)
    assert_cos_opt(opt, ref_cycle)

    # ap = animate(opt)
    # plt.show()


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


def test_spline_hei():
    kwargs = copy.copy(KWARGS)
    neb = NEB(get_geoms(), **NEB_KWARGS)
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
    dimer_result = dimer_method(geoms, **dimer_kwargs)
    dimer_cycles = dimer_result.dimer_cycles
    last_cycle = dimer_cycles[-1]
    ts_coords = last_cycle.trans_coords[1]
    print(f"Optimized TS coord: {ts_coords}")

    return opt



@pytest.mark.sd
def test_fix_first_neb():
    # First image is fixed at a non equilibrium geometry.
    coords = np.array(((-0.916, 1.034, 0), (1.94101, 3.85427, 0)))
    kwargs = copy.copy(KWARGS)
    neb = NEB(get_geoms(coords), fix_first=True, **NEB_KWARGS)
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 29)  # k = 0.01

    return opt


@pytest.mark.sd
def test_fix_last_neb():
    # Last image is fixed at a non equilibrium geometry.
    coords = np.array(((-1.05274, 1.02776, 0), (1.85, 3.57, 0)))
    kwargs = copy.copy(KWARGS)
    neb = NEB(get_geoms(coords), fix_last=True, **NEB_KWARGS)
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 28)  # k = 0.01

    return opt


@pytest.mark.sd
def test_fix_displaced_ends_neb():
    coords = np.array(((-0.916, 1.034, 0), (1.85, 3.57, 0)))
    kwargs = copy.copy(KWARGS)
    neb = NEB(get_geoms(coords), fix_ends=True, **NEB_KWARGS)
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 25)  # k = 0.01

    return opt


@pytest.mark.sd
def test_fix_end_immediate_climbing_neb():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 10
    kwargs["rms_force"] = 4e-2

    cos_kwargs = {
        "climb": True,
        "climb_rms": -1,
        "fix_ends": True,
    }
    neb = NEB(get_geoms(), **cos_kwargs, **NEB_KWARGS)
    opt = run_cos_opt(neb, SteepestDescent, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 39)

    return opt


@pytest.mark.fire
def test_fire_climb_neb():
    kwargs = copy.copy(KWARGS)
    kwargs["images"] = 10
    kwargs["dt_max"] = 0.2
    kwargs["climb"] = True
    kwargs["climb_rms"] = 0.0065
    neb = NEB(get_geoms(), fix_ends=True, **NEB_KWARGS)
    opt = run_cos_opt(neb, FIRE, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 38)  # k = 0.01

    return opt


@pytest.mark.skip
def test_fix_end_climbing_bfgs_neb():
    kwargs = copy.copy(KWARGS)
    # kwargs["bt_disable"] = True
    neb = NEB(get_geoms(), climb=True, fix_ends=True, **NEB_KWARGS)
    opt = run_cos_opt(neb, BFGS, **kwargs)

    assert(opt.is_converged)
    assert(opt.cur_cycle == 32)

    return opt


@pytest.mark.modified_broyden
def test_modified_broyden():
    pot = AnaPot()
    geom = pot.get_geom((-0.8333, 2, 0))
    force_getter = lambda x: geom.get_energy_and_forces_at(x)["forces"]
    def restrict_step(x, step, max_len=0.5):
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
        converged = norm < 1e-12
        if converged:
            print("Converged")
            break
    assert converged
    pot.plot()
    ax = pot.ax
    coords = np.array(coords)
    xs, ys = coords.T[:2]
    ax.plot(xs, ys, "ro-")
    # plt.show()
