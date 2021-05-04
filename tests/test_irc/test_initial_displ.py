import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.irc.initial_displ import cubic_displ_for_geom
from pysisyphus.irc import EulerPC


@pytest.fixture
def anapot_ts():
    ts_coords = (0.61173113, 1.49297317, 0.0)
    geom = AnaPot.get_geom(ts_coords)
    return geom


def test_cubic_displ(anapot_ts):
    step_plus, step_minus, _ = cubic_displ_for_geom(anapot_ts)
    ref_step_plus = (-0.013523096972, -0.009151256346, 0.0)
    np.testing.assert_allclose(step_plus, ref_step_plus)


def test_irc_cubic_displ(displ, anapot_ts):
    irc_kwargs = {
        "displ": "energy_cubic",
        "displ_energy": 1,
    }
    irc = EulerPC(anapot_ts, **irc_kwargs)
    irc.run()

    # calc = anapot_ts.calculator.plot_irc(irc, show=True)
