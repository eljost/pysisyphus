import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.MullerBrownSympyPot import MullerBrownPot
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.StringOptimizer import StringOptimizer
from pysisyphus.tsoptimizers.RSIRFOptimizer import RSIRFOptimizer


@pytest.mark.parametrize(
    "keep_last, ref_cycle", [
        (0, 9),
        (2, 11),
        (4, 13),
        (6, 15),
        (8, 17),
    ]
)
def test_anapot_growing_string(keep_last, ref_cycle):
    initial = AnaPot.get_geom((-1.05274, 1.02776, 0))
    final = AnaPot.get_geom((1.94101, 3.85427, 0))
    geoms = (initial, final)
    gs_kwargs = {
        "perp_thresh": 0.5,
        "reparam_check": "rms",
    }
    gs = GrowingString(geoms, lambda: AnaPot(), **gs_kwargs)

    opt_kwargs = {
        # "stop_in_when_full": 3,
        "stop_in_when_full": keep_last,
        "keep_last": keep_last,
    }
    opt = StringOptimizer(gs, **opt_kwargs)
    opt.run()

    # calc = AnaPot()
    # calc.anim_opt(opt, show=True)

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle


@pytest.mark.parametrize(
    "gs_kwargs_, opt_ref_cycle, tsopt_ref_cycle", [
        ({"climb": True, "climb_rms": 0.5, }, 22, 4),
        ({}, 22, 5),
    ]
)
def test_growing_string_climbing(gs_kwargs_, opt_ref_cycle, tsopt_ref_cycle):
    calc = AnaPot()
    geoms = calc.get_path(num=2)
    gs_kwargs = {
        "perp_thresh": 0.5,
        "reparam_check": "rms",
    }
    gs_kwargs.update(gs_kwargs_)
    cos = GrowingString(geoms, lambda: AnaPot(), **gs_kwargs)

    opt_kwargs = {
        "gamma_mult": True,
        "max_step": 0.04,
        "rms_force": 0.1,
        "rms_force_only": True,
    }
    opt = StringOptimizer(cos, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == opt_ref_cycle

    # calc = AnaPot()
    # calc.anim_opt(opt, show=True)

    hei_geom = cos.images[cos.get_hei_index()]
    tsopt = RSIRFOptimizer(hei_geom, thresh="gau_vtight")
    tsopt.run()

    assert tsopt.is_converged
    assert tsopt.cur_cycle == tsopt_ref_cycle


@pytest.mark.parametrize(
    "double_damp, ref_cycle", [
        pytest.param(False, 67, marks=pytest.mark.xfail),
        (True, 67),
    ]
)
def test_mullerbrown_string(double_damp, ref_cycle):
    calc = MullerBrownPot()
    geoms = calc.get_path(num=17)

    gs_kwargs = {
        "perp_thresh": 5.0,
        "max_nodes": 15,
    }
    # Reuse same calculator throughout, as the sympy call takes a while...
    gs = GrowingString(geoms, lambda: calc, **gs_kwargs)

    opt_kwargs = {
        "keep_last": 10,
        "lbfgs_when_full": True,
        "gamma_mult": True,
        "max_step": 0.03,
        "max_cycles": 350,
        "scale_step": "per_image",
        "double_damp": double_damp,
    }
    opt = StringOptimizer(gs, **opt_kwargs)
    opt.run()

    # calc.anim_opt(opt, show=True)

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle


@pytest.mark.parametrize(
    "double_damp, ref_cycle", [
        (False, 75),
        (True, 48),
    ]
)
def test_energy_reparametrization(double_damp, ref_cycle):
    calc = MullerBrownPot()
    geoms = calc.get_path(num=17)

    gs_kwargs = {
        "perp_thresh": 2.5,
        "max_nodes": 15,
        "param": "energy",
    }
    # Reuse same calculator throughout, as the sympy call takes a while...
    gs = GrowingString(geoms, lambda: calc, **gs_kwargs)

    opt_kwargs = {
        "keep_last": 10,
        "lbfgs_when_full": True,
        "max_step": 0.02,
        "max_cycles": 200,
        "rms_force": 7.5,
        "rms_force_only": True,
        "double_damp": double_damp
    }
    opt = StringOptimizer(gs, **opt_kwargs)
    opt.run()

    # calc.anim_opt(opt, show=True, energy_profile=True)

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle


def test_conjugate_gradient():
    calc = AnaPot()
    geoms = calc.get_path(2)
    gs_kwargs = {
        "perp_thresh": 0.5,
        "reparam_check": "rms",
    }
    gs = GrowingString(geoms, lambda: AnaPot(), **gs_kwargs)

    opt_kwargs = {
        "keep_last": 0,
        "rms_force": 0.02,
        "rms_force_only": True,
    }
    opt = StringOptimizer(gs, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == 23

    # calc.anim_opt(opt, show=True)


def test_climb_lanczos():
    calc = AnaPot()
    geoms = calc.get_path(2)
    gs_kwargs = {
        "perp_thresh": 0.5,
        "reparam_check": "rms",
        "climb": True,
        "climb_rms": 0.2,
        "climb_lanczos": True,
        "climb_lanczos_rms": 0.2,
    }
    gs = GrowingString(geoms, lambda: AnaPot(), **gs_kwargs)

    opt_kwargs = {
        "keep_last": 0,
        "rms_force": 0.02,
        "rms_force_only": True,
    }
    opt = StringOptimizer(gs, **opt_kwargs)
    opt.run()

    # calc.anim_opt(opt, show=True)

    assert opt.is_converged
    assert opt.cur_cycle == 23
