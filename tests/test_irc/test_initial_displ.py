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
    step = cubic_displ_for_geom(anapot_ts)
    print(np.array2string(step, precision=12))

    ref_step = (-0.013523096972, -0.009151256346, 0.0)
    np.testing.assert_allclose(step, ref_step)


@pytest.mark.parametrize(
    "displ", [
        "energy_cubic",
    ]
)
def test_irc_cubic_displ(displ, anapot_ts):
    irc_kwargs = {
        "displ": displ,
    }
    irc = EulerPC(anapot_ts, **irc_kwargs)
    irc.run()
