import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.cos.GrowingString import GrowingString
from pysisyphus.optimizers.StringOptimizer import StringOptimizer


@pytest.mark.parametrize(
    "keep_last, ref_cycle", [
        (0, 13),
        (2, 16),
        (4, 14),
        (6, 22),
        (8, 17),
        # (10, 17),
        # (12, 17),
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
        "stop_in_when_full": 3,
        "keep_last": keep_last,
    }
    opt = StringOptimizer(gs, **opt_kwargs)
    opt.run()

    assert opt.is_converged
    assert opt.cur_cycle == ref_cycle

    # calc = AnaPot()
    # calc.anim_opt(opt, show=True)
