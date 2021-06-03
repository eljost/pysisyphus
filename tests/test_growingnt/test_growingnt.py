import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.MullerBrownSympyPot import MullerBrownPot
from pysisyphus.cos.GrowingNT import GrowingNT
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS


def test_mb_growingnt():
    geoms = MullerBrownPot().get_minima()
    geoms = geoms[1], geoms[0]
    geom0, geom1 = geoms

    gnt_kwargs = {
        "step_len": 0.05,
        "final_geom": geom1,
        "rms_thresh": 0.5,
    }
    gnt = GrowingNT(geom0, **gnt_kwargs)
    opt_kwargs = {
        "max_cycles": 100,
    }
    opt = PreconLBFGS(gnt, **opt_kwargs)
    opt.run()

    calc = geoms[0].calculator
    calc.plot_geoms(gnt.images, show=True)
