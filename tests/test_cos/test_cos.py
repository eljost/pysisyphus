import matplotlib.pyplot as plt
import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.NFK import NFK
from pysisyphus.calculators.MullerBrownSympyPot import MullerBrownPot
from pysisyphus.cos.NEB import NEB
from pysisyphus.cos.SimpleZTS import SimpleZTS
from pysisyphus.Geometry import Geometry
from pysisyphus.interpolate.Interpolator import Interpolator
from pysisyphus.optimizers.LBFGS import LBFGS
from pysisyphus.optimizers.ConjugateGradient import ConjugateGradient
from pysisyphus.optimizers.QuickMin import QuickMin
from pysisyphus.optimizers.FIRE import FIRE
from pysisyphus.optimizers.SteepestDescent import SteepestDescent
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.plotters.AnimPlot import AnimPlot
from pysisyphus.run import run_from_dict
from pysisyphus.testing import using


def get_geoms():
    initial = AnaPot.get_geom((-1.05274, 1.02776, 0))
    final = AnaPot.get_geom((1.94101, 3.85427, 0))
    geoms = (initial, final)
    return geoms


def run_cos_opt(geoms, between, calc_cls, cos_cls, cos_kwargs, opt_cls, opt_kwargs):
    interpol = Interpolator(geoms, between=between)
    images = interpol.interpolate_all()

    for image in images:
        image.set_calculator(calc_cls())

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
        (SteepestDescent, {}, {}, 30, 5),
        (SteepestDescent, {}, {}, 32, 10),
        (ConjugateGradient, {}, {}, 44, 5),
        (QuickMin, {"dt": 0.1}, {}, 27, 5),
        (FIRE, {"dt_max": 0.2}, {}, 42, 5),
        (LBFGS, {"gamma_mult": True}, {}, 12, 5),
    ],
)
def test_anapot_neb(opt_cls, opt_kwargs_, neb_kwargs_, ref_cycle, between):
    geoms = get_geoms()

    neb_kwargs = {
        "k_min": 0.01,
    }
    neb_kwargs.update(neb_kwargs_)

    opt_kwargs = {}
    opt_kwargs.update(opt_kwargs_)

    opt = run_cos_opt(geoms, between, AnaPot, NEB, neb_kwargs, opt_cls, opt_kwargs)
    # ap = animate(opt)
    # plt.show()
    assert_cos_opt(opt, ref_cycle)


@pytest.mark.parametrize(
    "between, param, ref_cycle",
    [
        (5, "equal", 49),
        (10, "equal", 49),
        (5, "energy", 41),
        (10, "energy", 46),
    ],
)
def test_anapot_szts(between, param, ref_cycle):
    geoms = get_geoms()

    szts_kwargs = {
        "param": param,
    }

    opt_cls = SteepestDescent
    opt_kwargs = {
        "max_cycles": 100,
    }
    opt = run_cos_opt(
        geoms, between, AnaPot, SimpleZTS, szts_kwargs, SteepestDescent, opt_kwargs
    )
    # ap = animate(opt)
    # plt.show()
    assert_cos_opt(opt, ref_cycle)


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
    ap = AnimPlot(
        AnaPot(),
        opt,
        xlim=xlim,
        ylim=ylim,
        levels=levels,
        energy_profile=False,
        colorbar=False,
        figsize=(8, 6),
        save=False,
        title=False,
    )
    ap.animate()
    return ap


@using("pyscf")
def test_hcn_neb():
    run_dict = {
        "preopt": {
            "max_cycles": 3,
        },
        "interpol": {
            "type": "idpp",
            "between": 3,
        },
        "cos": {
            "type": "neb",
        },
        "opt": {
            "type": "qm",
            "align": True,
        },
        "calc": {
            "type": "pyscf",
            "pal": 2,
            "basis": "321g",
        },
        "geom": {
            "type": "cart",
            "fn": ["lib:hcn.xyz", "lib:hcn_iso_ts.xyz", "lib:nhc.xyz"],
        },
    }
    results = run_from_dict(run_dict)

    assert results.cos_opt.is_converged
    assert results.cos_opt.cur_cycle == 18


@pytest.mark.parametrize(
    "neb_kwargs, ref_cycle",
    [
        ({}, 34),
        ({"variable_springs": True, "k_min": 0.01, "k_max": 5}, 33),
    ],
)
def test_neb_springs(neb_kwargs, ref_cycle):
    calc = AnaPot()
    geoms = calc.get_path(15)
    neb = NEB(geoms, **neb_kwargs)
    opt = SteepestDescent(neb)
    opt.run()

    # calc.anim_opt(opt, show=True)

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle


@pytest.mark.parametrize(
    "k, ref_cycle",
    [
        (500, 70),
        (1000, 56),
        (2000, 62),
        (4000, 78),
    ],
)
def test_mullerbrown_neb(k, ref_cycle):
    geoms = MullerBrownPot().get_path(num=17)
    cos = NEB(geoms, k_max=k, k_min=k)

    opt_kwargs = {
        "max_step": 0.04,
        "gamma_mult": True,
        "keep_last": 10,
        "max_cycles": 100,
    }
    opt = LBFGS(cos, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle

    # calc = MullerBrownPot()
    # calc.anim_opt(opt, show=True)


def test_stiff_neb_nfk():
    nfk = NFK()
    # Coordinates from the paper
    # geoms = [nfk.get_geom(coords) for coords in ((-0.5, 0.25, 0.0), (2.7, -0.15, 0.0))]
    geoms = [nfk.get_geom(coords) for coords in ((-2.7, 0.13, 0.0), (2.7, -0.15, 0.0))]

    neb_kwargs = {
        "k_min": 10,
        "k_max": 10,
        # "bandwidth": 0.2,
    }

    opt_kwargs = {}
    between = 16
    # opt = run_cos_opt(geoms, between, NFK, NEB, neb_kwargs, LBFGS, opt_kwargs)
    opt = run_cos_opt(geoms, between, NFK, NEB, neb_kwargs, FIRE, opt_kwargs)

    # ap = AnimPlot(nfk, opt, xlim=nfk.xlim, ylim=nfk.ylim)
    # ap.animate()
    # plt.show()


@pytest.mark.parametrize("climb, ref_cycles", ((True, 209), ("one", 109), (False, 120)))
def test_neb_climb(climb, ref_cycles):
    geoms = get_geoms()

    neb_kwargs = {
        "variable_springs": True,
        "climb": climb,
    }

    between = 3
    opt_kwargs = {
        "max_cycles": 500,
        "thresh": "gau",
        "max_step": 0.05,
    }
    opt_cls = ConjugateGradient
    opt = run_cos_opt(geoms, between, AnaPot, NEB, neb_kwargs, opt_cls, opt_kwargs)
    # ap = animate(opt)
    # plt.show()
    assert opt.is_converged
    assert opt.cur_cycle == ref_cycles
