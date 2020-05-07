#!/usr/bin/env python3

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.cos.AdaptiveNEB import AdaptiveNEB
from pysisyphus.optimizers.ConjugateGradient import ConjugateGradient


def test_anapot_aneb():
    image_num = 10
    calc = AnaPot()
    all_geoms = calc.get_path(image_num)
    aneb_kwargs = {
        # "keep_hei": True,
    }
    aneb = AdaptiveNEB(all_geoms, **aneb_kwargs)

    opt_kwargs = {
        "rms_force": 3e-3,
    }
    opt = ConjugateGradient(aneb, **opt_kwargs)
    opt.run()

    ap = calc.anim_opt(opt, show=True)

    assert opt.is_converged
    assert opt.cur_cycle == 20
