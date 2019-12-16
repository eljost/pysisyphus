#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.irc import *


def plot_irc(irc, title=None):
    geom = irc.geometry
    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    ax.plot(*irc.all_coords.T[:2], "ro-")
    if title:
        ax.set_title(title)
    plt.show()


@pytest.mark.parametrize(
    "irc_cls, mod_kwargs, ref", [
        (DampedVelocityVerlet, {}, None),
        # (Euler, {}, None),
        # (EulerPC, {}, None),
        # (GonzalesSchlegel, {}, None),
        # (IMKMod, {}, None),
        # (RK4, {}, None),
        # (LQA, {}, None),
    ]
)
def test_anapot_irc(irc_cls, mod_kwargs, ref):
    geom = AnaPot().get_geom((0.61173, 1.49297, 0.))

    kwargs = {
        "step_length": 0.05,
        "rms_grad_thresh": 1e-2,
    }
    kwargs.update(**mod_kwargs)

    irc = irc_cls(geom, **kwargs)
    irc.run()

    fc = irc.all_coords[0]
    bc = irc.all_coords[-1]
    forward_ref = np.array((-1.0527, 1.0278,  0.))
    backward_ref = np.array((1.941, 3.8543, 0.))
    assert np.allclose(fc, forward_ref, atol=0.03)
    assert np.allclose(bc, backward_ref, atol=0.03)

    plot_irc(irc, irc.__class__.__name__)
