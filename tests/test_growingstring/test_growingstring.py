import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.MullerBrownSympyPot import MullerBrownPot
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.StringOptimizer import StringOptimizer


@pytest.mark.parametrize(
    "keep_last, ref_cycle", [
        (0, 10),
        (2, 12),
        (4, 14),
        (6, 16),
        (8, 18),
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
        "gamma": 10.,
        # "stop_in_when_full": 3,
        "stop_in_when_full": keep_last,
        "keep_last": keep_last,
    }
    opt = StringOptimizer(gs, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle

    # calc = AnaPot()
    # calc.anim_opt(opt, show=True)


def test_mullerbrown_string():
    calc = MullerBrownPot()
    geoms = calc.get_path(num=17)

    gs_kwargs = {
        "perp_thresh": 5.0,
        "reparam_check": "rms",
        "max_nodes": 15,
    }
    # Reuse same calculator throughout, as the sympy call takes a while...
    gs = GrowingString(geoms, lambda: calc, **gs_kwargs)

    opt_kwargs = {
        "keep_last": 10,
        "lbfgs_when_full": True,
        "gamma_mult": True,
        "max_step": 0.04,
        "max_cycles": 75,
    }
    opt = StringOptimizer(gs, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == 67

    # calc = MullerBrownPot()
    # calc.anim_opt(opt, show=True)


def test_gs():
    calc = AnaPot()
    geoms = calc.get_path(num=2)
    gs_kwargs = {
        "perp_thresh": 0.5,
        "reparam_check": "rms",
    }
    # gs_kwargs = {}
    gs = GrowingString(geoms, lambda: AnaPot(), **gs_kwargs)

    opt_kwargs = {
        "gamma": 10.,
        "gamma_mult": True,
        "max_step": 0.04,
        "rms_force": 0.1,
        "rms_force_only": True,
    }
    opt = StringOptimizer(gs, **opt_kwargs)
    opt.run()

    # assert opt.is_converged
    # assert opt.cur_cycle == ref_cycle

    # calc = AnaPot()
    # calc.anim_opt(opt, show=True)
