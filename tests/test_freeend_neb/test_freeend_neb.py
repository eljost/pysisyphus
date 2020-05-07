#!/usr/bin/env python3

import copy

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pysisyphus.plotters.AnimPlot import AnimPlot
from pysisyphus.calculators.FreeEndNEBPot import FreeEndNEBPot
from pysisyphus.cos.FreeEndNEB import FreeEndNEB
from pysisyphus.cos.NEB import NEB
from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.SteepestDescent import SteepestDescent
from pysisyphus.optimizers.QuickMin import QuickMin
# from pysisyphus.optimizers.BFGS import BFGS

KWARGS = {
    "images": 21,
    "max_cycles": 100,
    "convergence": {
        "max_force_thresh": 6e-3,
        "rms_force_thresh": 5e-3,
        "max_step_thresh": 6e-4,
        "rms_step_thresh": 3e-4,
    },
    "dump": False,
}


@pytest.mark.skip
def test_steepest_descent_neb():
    kwargs = copy.copy(KWARGS)
    calc = FreeEndNEBPot()
    all_geoms = calc.get_path(50)[:15]
    # from pysisyphus.cos.NEB import NEB
    neb_kwargs = {
        "fix_ends": False,
        "k_min": 50,
        "k_max": 53,
        "mod": True,
    }
    neb = FreeEndNEB(all_geoms, **neb_kwargs)
    # neb = NEB(all_geoms)#, k_min=50, k_max=53)#, mod=True)#, k_min=50, k_max=50.3)
    # opt = QuickMin(neb, max_cycles=60)
    opt = SteepestDescent(neb, max_cycles=3)
    opt.run()
    calc.anim_opt(opt, show=True)

    return opt


@pytest.mark.skip
def yolo():
    from pysisyphus.calculators.LEPSBase import LEPSBase
    lb = LEPSBase(pot_type="harmonic")
    initial = np.array((0.7423, 1.18342, 0))
    final = np.array((2.7234, -1.08509, 0))
    coords = (initial, final)
    atoms = ("H")
    geoms = [Geometry(atoms, c) for c in coords]
    cos = NEB(geoms)
    cos.interpolate(21)
    opt = SteepestDescent(cos, max_cycles=10)
    for img in cos.images:
        img.set_calculator(LEPSBase("harmonic"))
    opt.run()
    levels = (-20, 20, 250)
    ap = AnimPlot(lb, opt, xlim=lb.xlim, ylim=lb.ylim, levels=levels, tight_layout=False)
    ap.animate()
    # plt.show()


@pytest.mark.skip
def test_mb():
    from pysisyphus.calculators.MullerBrownPot import MullerBrownPot
    initial = np.array((-0.5588, 1.437, 0))
    final = np.array((-0.22404, 1.0055, 0))
    calc = MullerBrownPot()
    geoms = [Geometry(("H", ), coords) for coords in (initial, final) ]
    from pysisyphus.interpolate import interpolate
    geoms = interpolate(*geoms, between=5)
    for geom in geoms:
        geom.set_calculator(calc)
    cos = FreeEndNEB(geoms, k_min=10, k_max=20)
    # cos = NEB(geoms)
    opt = QuickMin(cos, max_cycles=100)
    opt.run()
    levels=(-150, -15, 40)
    xlim = (-1, 0.2)
    ylim = (0.4, 1.8)
    mb = MullerBrownPot()
    ap = AnimPlot(mb, opt, xlim=xlim, ylim=ylim, levels=levels)
    ap.animate()
    # plt.show()


if __name__ == "__main__":
    # opt = test_steepest_descent_neb()
    test_mb()
    # yolo()
