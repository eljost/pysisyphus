#!/usr/bin/env python3

import numpy as np

from pysisyphus.calculators.MullerBrownSympyPot2D import MullerBrownSympyPot2D
from pysisyphus.Geometry import Geometry
from pysisyphus.irc.DampedVelocityVerlet import DampedVelocityVerlet
from pysisyphus.irc.Euler import Euler
from pysisyphus.irc.GonzalesSchlegel import GonzalesSchlegel
from pysisyphus.irc.PlotAnaPot import PlotAnaPot


def prepare():
    atoms = ("H", )
    # TS coordinates
    calc, ts_coords = (MullerBrownSympyPot2D(), np.array((-0.822, 0.624)))
    geometry = Geometry(atoms, ts_coords)
    geometry.set_calculator(calc)

    xlim = (-1.75, 1.25)
    ylim = (-0.5, 2.25)
    levels=(-150, -15, 80)
    plotter = PlotAnaPot(geometry, xlim, ylim, levels)

    return geometry, plotter


def test_mullerbrown_gs_irc():
    geometry, plotter = prepare()

    gs = GonzalesSchlegel(geometry, step_length=0.3, max_steps=6)
    gs.run()
    assert(gs.forward_step == 5)
    assert(gs.backward_step == 4)

    plotter.plot_gs(gs.all_coords, gs.pivot_coords, gs.micro_coords)
    #plotter.show()


def test_mullerbrown_dvv_irc():
    geometry, plotter = prepare()
    dvv = DampedVelocityVerlet(geometry, v0=0.04, max_steps=71)
    dvv.run()
    assert(dvv.forward_step == 36)
    assert(dvv.backward_step == 70)

    plotter.plot(dvv.all_coords)
    #plotter.show()


def test_mullerbrown_euler_irc():
    geometry, plotter = prepare()
    euler = Euler(geometry, max_steps=150, step_length=0.05)
    euler.run()
    assert(euler.forward_step == 22)
    assert(euler.backward_step == 17)

    plotter.plot(euler.all_coords)
    #plotter.show()


def test_mullerbrown_euler_no_mw_irc():
    geometry, plotter = prepare()
    euler = Euler(geometry, max_steps=150, step_length=0.05, mass_weight=False)
    euler.run()
    #assert(euler.forward_step == 22)
    #assert(euler.backward_step == 17)

    plotter.plot(euler.all_coords)
    #plotter.show()


if __name__ == "__main__":
    test_mullerbrown_gs_irc()
    #test_mullerbrown_dvv_irc()
    #test_mullerbrown_euler_irc()
    #test_mullerbrown_euler_no_mw_irc()
