#!/usr/bin/env python3

import itertools as it

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.AnaPotCBM import AnaPotCBM
from pysisyphus.calculators.Gaussian16 import Gaussian16
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_library
from pysisyphus.tsoptimizers.dimer import dimer_method


def get_geoms(coords=None):
    if coords is None:
        # left = np.array((0.188646, 1.45698, 0))
        # right = np.array((0.950829, 1.54153, 0))
        left = np.array((0.354902, 1.34229, 0))
        right = np.array((0.881002, 1.71074, 0))
        right = np.array((0.77, 1.97, 0))

        # Very close
        # left = np.array((0.531642, 1.41899, 0))
        # right = np.array((0.702108, 1.57077, 0))
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


def plot_dimer_cycles(dimer_cycles, pot, true_ts=None):
    pot.plot()

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
    if true_ts:
        ts_x, ts_y = true_ts
        ax.scatter(ts_x, ts_y, s=20, zorder=25, c="k")
    pot.ax.legend()
    plt.show()


def test_anapot(trans_opt, trans_memory):
    geoms = get_geoms()
    calc_getter = AnaPot
    dimer_kwargs = {
        "ana_2dpot": True,
        "restrict_step": "max",
        "angle_tol": 0.5,
        "f_thresh": 1e-4,
        "rot_opt": "mb",
        "trans_opt": trans_opt,
        "trans_memory": trans_memory,
        "f_tran_mod": False,
    }
    result = dimer_method(geoms, calc_getter, **dimer_kwargs)
    return result


def anapot_tester():
    trans_opts = ("lbfgs", "mb")
    trans_memories = range(2, 8)
    # trans_opts = ("mb", )
    # trans_memories = (4, )
    results = dict()
    true_ts = (0.61173, 1.49297, 0.)
    ref_cycles = {
        'lbfgs_2': 9,
        'lbfgs_3': 9,
        'lbfgs_4': 7,
        'lbfgs_5': 7,
        'lbfgs_6': 7,
        'lbfgs_7': 7,
        'mb_2': 5,
        'mb_3': 5,
        'mb_4': 5,
        'mb_5': 5,
        'mb_6': 5,
        'mb_7': 5,
    }
    for to, tm in it.product(trans_opts, trans_memories):
        dimer_result = test_anapot(to, tm)
        key = f"{to}_{tm}"
        results[key] = len(dimer_result.dimer_cycles)
        np.testing.assert_allclose(dimer_result.geom0.coords, true_ts, atol=1e-4)
        assert ref_cycles[key] == len(dimer_result.dimer_cycles)
    print(results)
    # plot_dimer_cycles(dimer_result.dimer_cycles, pot=AnaPot(), true_ts=true_ts[:2])



def plot_anapotcbm_curvature():
    pot = AnaPotCBM()
    pot.plot()
    xs = np.linspace(-1.25, 1.25, num=50)
    # ys = np.linspace(-0.75, 0.75, num=50)
    ys = np.linspace(-1, 1, num=50)
    X, Y = np.meshgrid(xs, ys)
    z = list()
    for x_, y_ in zip(X.flatten(), Y.flatten()):
        g = pot.get_geom((x_, y_, 0))
        H = g.hessian
        w, v = np.linalg.eigh(H)
        z.append(
            1 if (w < 0).any() else 0
        )
    Z = np.array(z).reshape(X.shape)
    ax = pot.ax
    ax.contourf(X, Y, Z, cmap=cm.Reds)#, alpha=0.5)
    plt.show()


def test_anapotcbm():
    calc_getter = AnaPotCBM
    # geom = AnaPotCBM().get_geom((0.818, 0.2233, 0.0))
    geom = AnaPotCBM().get_geom((0.2, 0.2, 0.0))
    # geom = AnaPotCBM().get_geom((0.5, 0.2, 0.0))
    geom = AnaPotCBM().get_geom((0.9, 0.8, 0.0))
    v, w = np.linalg.eigh(geom.hessian)
    N_imag = w[:,0]
    geoms = [geom, ]
    dimer_kwargs = {
        "ana_2dpot": True,
        "restrict_step": "max",
        "N_init": N_imag,
        "trans_opt": "mb",
    }
    true_ts = (0, 0)
    dimer_result = dimer_method(geoms, calc_getter, **dimer_kwargs)
    dimer_cycles = dimer_result.dimer_cycles
    plot_dimer_cycles(dimer_cycles, pot=AnaPotCBM(), true_ts=true_ts)


def test_hcn_iso_dimer(trans_opt, trans_memory):

    calc_kwargs = {
        "route": "PM6",
        "pal": 4,
        "mem": 1000,
    }
    def calc_getter():
        return Gaussian16(**calc_kwargs)

    geom = geom_from_library("hcn_iso_pm6_near_ts.xyz")
    geom.set_calculator(calc_getter())
    geoms = [geom, ]

    N_init = np.array(
        (0.6333, 0.1061, 0.5678, 0.171, 0.11, 0.3373, 0.0308, 0.1721, 0.282)
    )
    dimer_kwargs = {
        "max_step": 0.04,
        "dR_base": 0.01,
        "N_init": N_init,
        # "rot_opt": "mb",
        "trans_opt": trans_opt,
        "trans_memory": trans_memory,
        "angle_tol": 5,
        "f_thresh": 1e-4,
    }
    dimer_result = dimer_method(geoms, calc_getter, **dimer_kwargs)
    return dimer_result


def hcn_tester():
    trans_opts = ("lbfgs", "mb")
    # trans_opts = ("mb", )
    trans_memories = range(3, 6)

    # Original f_trans
    ref_evals = {
        'lbfgs_3': 18,
        'lbfgs_4': 18,
        'lbfgs_5': 18,
        'mb_3': 16,
        'mb_4': 16,
        'mb_5': 16,
    }
    # Modified f_trans
    ref_evals = {
        'lbfgs_3': 18,
        'lbfgs_4': 18,
        'lbfgs_5': 18,
        'mb_3': 16,
        'mb_4': 18,
        'mb_5': 18,
    }
    results = dict()
    for to, tm in it.product(trans_opts, trans_memories):
        dimer_result = test_hcn_iso_dimer(to, tm)
        key = f"{to}_{tm}"
        results[key] = dimer_result.force_evals
        assert ref_evals[key] == dimer_result.force_evals
    print(results)


if __name__ == "__main__":
    anapot_tester()
    hcn_tester()
    # plot_anapotcbm_curvature()
    # test_anapotcbm()
