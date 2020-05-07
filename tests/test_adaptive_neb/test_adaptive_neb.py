#!/usr/bin/env python3

import copy

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pysisyphus.plotters.AnimPlot import AnimPlot
from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.cos.AdaptiveNEB import AdaptiveNEB
from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.SteepestDescent import SteepestDescent


KWARGS = {
    "images": 3,
    "max_cycles": 50,
    "dump": False,
}


def run_cos_opt(cos, Opt, images, **kwargs):
    opt = Opt(cos, **kwargs)
    for img in cos.images:
        img.set_calculator(AnaPot())
    opt.run()

    return opt


def animate_bare(opt):
    xlim = (-2, 2.5)
    ylim = (0, 5)
    levels = (-3, 4, 80)
    ap = AnimPlot(AnaPot(), opt, xlim=xlim, ylim=ylim, levels=levels,
                  energy_profile=False, colorbar=False, figsize=(8, 6))
    ap.animate()
    return ap


@pytest.mark.sd
def test_steepest_descent_neb():
    kwargs = copy.copy(KWARGS)
    all_geoms = AnaPot().get_path(5)
    aneb = AdaptiveNEB(all_geoms)#, keep_hei=False)
    opt = run_cos_opt(aneb, SteepestDescent, **kwargs)

    return opt


if __name__ == "__main__":
    opt = test_steepest_descent_neb()
    ap = animate_bare(opt)
    plt.show()
