#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.Geometry import Geometry
from pysisyphus.tsoptimizers.dimer import dimer_method

def get_geoms(coords=None):
    if coords is None:
        # left = np.array((0.188646, 1.45698, 0))
        # right = np.array((0.950829, 1.54153, 0))
        # left = np.array((0.354902, 1.34229, 0))
        # right = np.array((0.881002, 1.71074, 0))
        left = np.array((0.531642, 1.41899, 0))
        right = np.array((0.702108, 1.57077, 0))
        coords = (right, left)

        # near_ts = np.array((0.553726, 1.45458, 0))
        # coords = (near_ts, )

        # left_far = np.array((-0.455116, 0.926978, 0))
        # right_far = np.array((-0.185653, 1.02486, 0))
        # coords = (left_far, right_far)

    atoms = ("H")
    geoms = [Geometry(atoms, c) for c in coords]
    for geom in geoms:
        geom.set_calculator(AnaPot())
    return geoms


def plot_dimer(dimer, ax, label=None, color=None, marker="o"):
    lines = ax.plot(dimer[:,0], dimer[:,1], marker=marker,
                    label=label, color=color)
    return lines


def plot_dimer_cycles(dimer_cycles):
    pot = AnaPot()
    levels = np.linspace(-2.8, 3.6, 50)
    pot.plot(levels=levels)

    ax = pot.ax
    for i, dc in enumerate(dimer_cycles):
        label = f"Cycle {i}"
        org_lbl = f"Org {i}"
        trial_lbl = f"Trial {i}"
        rot_lbl = f"Rot {i}"
        # org = plot_dimer(dc.org_coords, ax, label=org_lbl)
        # color = org[0].get_color()
        # trial = plot_dimer(dc.trial_coords, ax,
                         # label=trial_lbl, color=color, marker="x")
        # rot = plot_dimer(dc.rot_coords, ax,
                         # label=rot_lbl, color=color, marker=".")
        rot = plot_dimer(dc.rot_coords, ax,
                         label=rot_lbl, marker=".")
    pot.ax.legend()
    plt.show()


def run():
    geoms = get_geoms()
    calc_getter = AnaPot
    dimer_kwargs = {
        "max_step": 0.1,
        "ana_2dpot": True,
    }
    dimer_cycles = dimer_method(geoms, calc_getter, **dimer_kwargs)
    plot_dimer_cycles(dimer_cycles[-5:])


def test_hcn_iso_dimer():
    from pysisyphus.helpers import geom_from_library
    from pysisyphus.calculators.Gaussian16 import Gaussian16

    calc_kwargs = {
        "route": "PM6",
        "pal": 4,
        "mem": 1000,
    }
    def calc_getter():
        return Gaussian16(**calc_kwargs)

    geom = geom_from_library("hcn_iso_ts_guess.xyz")
    geom = geom_from_library("hcn_iso_pm6_ts.xyz")
    geom = geom_from_library("hcn_iso_pm6_near_ts.xyz")
    geom.set_calculator(calc_getter())
    geoms = [geom, ]

    dimer_kwargs = {
        #"max_step": 0.1,
        "ana_2dpot": False,
    }
    dimer_cycles = dimer_method(geoms, calc_getter, **dimer_kwargs)
    # plot_dimer_cycles(dimer_cycles[-5:])


if __name__ == "__main__":
    run()
    # test_hcn_iso_dimer()
