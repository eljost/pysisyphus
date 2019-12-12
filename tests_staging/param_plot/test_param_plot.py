import os
import pathlib

import numpy as np

from pysisyphus.irc.ParamPlot import ParamPlot

THIS_DIR = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))

def run_param_plot(coords, p1inds, p2inds, prefix):
    param_plot = ParamPlot(coords, p1inds, p2inds)
    param_plot.plot()
    #param_plot.show()
    param_plot.save(THIS_DIR, prefix)

    return param_plot

def test_hcn_iso_hfsto3g():
    coords_fn = THIS_DIR / "hcn_iso_hfsto3g_gs_0_2_.coords"
    coords = np.loadtxt(coords_fn)
    prefix = "hcn_iso_hfsto3g_gs_0_2"
    p1inds = (0, 1)
    p2inds = (1, 0, 2)
    param_plot = run_param_plot(coords, p1inds, p2inds, prefix)
