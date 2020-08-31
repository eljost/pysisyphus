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


@pytest.mark.skip
def test_anapot_growing_string():
    initial = MullerBrownPot.get_geom((-0.5592, 1.443, 0))
    final = MullerBrownPot.get_geom((0.605, 0.036, 0))
    calc = MullerBrownPot()
    geoms = (initial, final)
    from pysisyphus.optimizers.RFOptimizer import RFOptimizer
    for geom in geoms:
        opt = RFOptimizer(geom, thresh="gau_tight", hessian_recalc=1,)
        opt.run()

    from pysisyphus.cos.NEB import NEB
    from pysisyphus.interpolate import interpolate
    geoms = interpolate(*geoms, 18)
    for geom in geoms:
        geom.set_calculator(calc)
    cos = NEB(geoms, k_max=5.1, k_min=5.)
    from pysisyphus.optimizers.QuickMin import QuickMin
    from pysisyphus.optimizers.LBFGS import LBFGS
    opt = QuickMin(cos, max_cycles=10)
    opt.run()
    opt = LBFGS(cos, max_step=0.02, gamma_mult=True)
    opt.run()
    calc = MullerBrownPot()
    calc.anim_opt(opt, show=True)
    import sys; sys.exit()

    gs_kwargs = {
        "perp_thresh": 25.0,
        "reparam_check": "rms",
        "max_nodes": 16,
        # "reparam_every": 10,
    }
    gs = GrowingString(geoms, lambda: MullerBrownPot(), **gs_kwargs)

    opt_kwargs = {
        # "gamma": 20,
        "gamma": 100,
        "lbfgs_when_full": False,
        "max_step": 0.03,
        "max_cycles": 200,
    }
    opt = StringOptimizer(gs, **opt_kwargs)
    opt.run()

    # import pdb; pdb.set_trace()
    calc = MullerBrownPot()
    calc.anim_opt(opt, show=True)
