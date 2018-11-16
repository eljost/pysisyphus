#!/usr/bin/env python3

import copy

import matplotlib.pyplot as plt
import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.MullerBrownSympyPot import MullerBrownPot
from pysisyphus.calculators.FourWellAnaPot import FourWellAnaPot
from pysisyphus.Geometry import Geometry
from pysisyphus.cos.GrowingNT import GrowingNT


def get_geoms(coords, calc_getter):
    atoms = ("H")
    geoms = [Geometry(atoms, c) for c in coords]
    for geom in geoms:
        geom.set_calculator(calc_getter())
    return geoms


def plot(gnt, calc, levels=None):
    calc.plot(levels)
    points = np.array(gnt.points)
    conv = np.array(gnt.conv_points)
    # for i, p in enumerate(points):
        # print(i, p)
    px = points[:,0]
    py = points[:,1]
    cx = conv[:,0]
    cy = conv[:,1]
    calc.ax.plot(px, py, "o-", c="r")
    calc.ax.plot(cx, cy, "X-", ms="5", c="k")
    plt.show()


def test_anapot_growingnt():
    coords = (
        (-1.05274, 1.02776, 0),
        (1.94101, 3.85427, 0),
    )
    levels = np.linspace(-3, 4, 80)
    calc_getter = AnaPot
    eps = .05
    damp = .05
    gnt = GrowingNT(get_geoms(coords, calc_getter), calc_getter,
                    eps=eps, damp=damp)

    gnt.run()
    plot(gnt, calc_getter(), levels)


def test_mullerbrown_growingnt():
    coords = (
        (0.614, 0.031, 0),
        (-.563, 1.43, 0),
    )
    levels=np.linspace(-150, -15, 40)
    eps = .008
    damp = .00065

    calc_getter = MullerBrownPot
    gnt = GrowingNT(get_geoms(coords, calc_getter), calc_getter, 
                    eps=eps, damp=damp, max_nodes=23)
    gnt.run()
    plot(gnt, calc_getter(), levels)


def test_four_well_growingnt():
    coords = (
        (1.124, -1.485, 0.0),
        (-1.174, 1.477, 0.0),
    )
    eps = .25
    damp = .05

    # See 10.1063/1.1885467 Sec. V.B., Eq. 11 and Fig 4a)
    calc_getter = FourWellAnaPot
    gnt = GrowingNT(get_geoms(coords, calc_getter), calc_getter, 
                    eps=eps, damp=damp, max_nodes=22, readjust=False)
    gnt.run()
    plot(gnt, calc_getter())

if __name__ == "__main__":
    # test_anapot_growingnt()
    # test_mullerbrown_growingnt()
    test_four_well_growingnt()
