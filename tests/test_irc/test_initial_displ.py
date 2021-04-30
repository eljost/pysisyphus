import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.irc.initial_displ import cubic_displ_for_geom


def test_cubic_displ():
    ts_coords = (0.61173113, 1.49297317, 0.0)
    geom = AnaPot.get_geom(ts_coords)
    step = cubic_displ_for_geom(geom)
    print(np.array2string(step, precision=12))

    ref_step = (-0.013523096972, -0.009151256346, 0.0)
    np.testing.assert_allclose(step, ref_step)
