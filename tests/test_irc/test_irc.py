#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.irc.Euler import Euler
from pysisyphus.irc.EulerPC import EulerPC
from pysisyphus.Geometry import Geometry
from pysisyphus.irc.GonzalesSchlegel import GonzalesSchlegel
from pysisyphus.irc.IMKMod import IMKMod
from pysisyphus.irc import RK4
# from pysisyphus.calculators.MullerBrownSympyPot import MullerBrownPot
# from pysi
from pysisyphus.calculators.MullerBrownSympyPot import MullerBrownPot


def get_geom():
    ts_coords = (0.61173, 1.49297, 0.)
    return AnaPot().get_geom(ts_coords)


def test_gs():
    geom = get_geom()
    irc_kwargs = {
        "step_length": 0.2,
        "max_steps": 50,
        "rms_grad_thresh": 1e-2,
    }
    irc = GonzalesSchlegel(geom, **irc_kwargs)
    irc.run()

    return irc


def test_imk():
    geom = get_geom()
    # irc_kwargs = {
        # "step_length": 0.4,
        # "max_steps": 50,
        # "rms_grad_thresh": 1e-2,
    # }
    irc_kwargs = {}
    irc = IMKMod(geom, **irc_kwargs)
    irc.run()

    return irc


def test_rk4():
    geom = AnaPot.get_geom((0.61173, 1.49297, 0))
    irc_kwargs = {
        "step_length": 0.2,
    }
    irc = RK4(geom, **irc_kwargs)
    irc.run()

    return irc


def plot_anapot_irc(irc):
    calc = AnaPot()
    calc.plot()
    fig, ax = calc.fig, calc.ax

    ax.scatter(0.61173, 1.49297, s=75, c="k")
    ax.scatter(-1.05274, 1.02776, s=75, c="k")
    ax.scatter(1.94101, 3.85427, s=75, c="k")

    def label_steps(ax, xs, ys):
        for i, (x, y) in enumerate(zip(xs, ys)):
            ax.annotate(f"{i}", (x, y))

    try:
        fw_coords = np.array(irc.forward_coords)
        ax.plot(fw_coords[:, 0], fw_coords[:, 1], "ro", ls="-", label="forward")
        label_steps(ax, fw_coords[:,0][::-1], fw_coords[:,1][::-1])
    except AttributeError:
        pass
    try:
        bw_coords = np.array(irc.backward_coords)
        ax.plot(bw_coords[:, 0], bw_coords[:, 1], "bo", ls="-", label="backward")
        label_steps(ax, bw_coords[:,0], bw_coords[:,1])
    except AttributeError:
        pass
    ax.legend()
    plt.show()


def test_lqa():
    from pysisyphus.irc.LQA import LQA
    ts_coords = (0.61173, 1.49297, 0.)
    geom = AnaPot.get_geom(ts_coords)
    irc = LQA(geom, step_length=.2)
    irc.run()

    coords = irc.all_coords
    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    ax.plot(*coords.T[:2], "ro-")
    plt.show()


def test_eulerpc():
    # ts_coords = (0.61173, 1.49297, 0.)
    # geom = AnaPot.get_geom(ts_coords)
    ts_coords = (-0.822, 0.624, 0.)
    geom = MullerBrownPot .get_geom(ts_coords)

    irc_kwargs = {
        # "step_length": 0.5,
        # "step_length": 1.5,
        # "step_length": 150,
        # "step_length": 1.5,
        # "step_length": .2,
        # "step_length": .4,
        # "step_length": .4,
        # "step_length": .3,
        "step_length": .1,
        "displ": "length",
        "displ_length": .05,
        # "hessian_update": "bfgs",
        "hessian_update": "bofill",
    }
    irc = EulerPC(geom, **irc_kwargs)
    irc.run()

    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    ax.plot(*irc.all_coords_umw.T[:2], "ro-")
    ax.set_xlim(-1.1, 0.1)
    ax.set_ylim( 0.3, 1.6)
    plt.show()


if __name__ == "__main__":
    # irc = test_imk()
    # irc = test_rk4()
    # plot_anapot_irc(irc)
    # irc = test_gs()
    # plot_anapot_irc(irc)
    # test_lqa()
    test_eulerpc()
