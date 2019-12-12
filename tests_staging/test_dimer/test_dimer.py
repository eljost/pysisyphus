#!/usr/bin/env python3

import itertools as it
from pprint import pprint

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.AnaPotCBM import AnaPotCBM
from pysisyphus.calculators.Gaussian16 import Gaussian16
from pysisyphus.calculators.XTB import XTB
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
        # "multiple_translations": True,
    }
    result = dimer_method(geoms, calc_getter, **dimer_kwargs)
    return result


def test_anapot_rot():
    geoms = get_geoms()
    calc_getter = AnaPot
    dimer_kwargs = {
        "ana_2dpot": True,
        "restrict_step": "max",
        "angle_tol": 0.5,
        "f_thresh": 1e-4,
        "rot_opt": "mb",
        "f_tran_mod": False,
        "rot_type": "direct",
        # "rot_f_thresh": 0.1,
        "multiple_translations": False,
        # "max_cycles": 1,
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
    """
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
    """
    ref_force_evals = {
        'lbfgs_2': 32,
        'lbfgs_3': 32,
        'lbfgs_4': 26,
        'lbfgs_5': 26,
        'lbfgs_6': 26,
        'lbfgs_7': 26,
        'mb_2': 21,
        'mb_3': 21,
        'mb_4': 21,
        'mb_5': 21,
        'mb_6': 21,
        'mb_7': 21
    }
    result = dict()
    for to, tm in it.product(trans_opts, trans_memories):
        dimer_result = test_anapot(to, tm)
        key = f"{to}_{tm}"
        result[key] = dimer_result.force_evals
        np.testing.assert_allclose(dimer_result.geom0.coords, true_ts, atol=1e-4)
        assert dimer_result.force_evals == ref_force_evals[key]
    # plot_dimer_cycles(dimer_result.dimer_cycles, pot=AnaPot(), true_ts=true_ts[:2])
    pprint(result)


def test_anapot_cbm_rot():
    pot = AnaPot()
    geoms = (pot.get_geom((0.44, 1.54, 0)), )
    N_init = (-0.2, 1, 0)
    calc_getter = AnaPot
    dimer_kwargs = {
        "ana_2dpot": True,
        "restrict_step": "max",
        "angle_tol": 0.5,
        "f_thresh": 1e-4,
        "rot_opt": "lbfgs",
        "trans_opt": "lbfgs",
        "trans_memory": 5,
        "f_tran_mod": False,
        "N_init": N_init,
        "rot_f_thresh": 1e-2,
        # "multiple_translations": True,
    }
    result = dimer_method(geoms, calc_getter, **dimer_kwargs)
    true_ts = (0.61173, 1.49297, 0.)
    # plot_dimer_cycles(result.dimer_cycles, pot=AnaPot(), true_ts=true_ts[:2])
    return result


def plot_anapotcbm_curvature():
    pot = AnaPotCBM()
    pot.plot()
    num = 25
    xs = np.linspace(-1.25, 1.25, num=num)
    # ys = np.linspace(-0.75, 0.75, num=50)
    ys = np.linspace(-1, 1, num=num)
    X, Y = np.meshgrid(xs, ys)
    z = list()
    neg = list()
    for x_, y_ in zip(X.flatten(), Y.flatten()):
        g = pot.get_geom((x_, y_, 0))
        H = g.hessian
        w, v = np.linalg.eigh(H)
        z.append(
            1 if (w < 0).any() else 0
        )
        if (w < 0).any():
            neg.append((x_, y_))
    Z = np.array(z).reshape(X.shape)
    ax = pot.ax
    # ax.contourf(X, Y, Z, cmap=cm.Reds)#, alpha=0.5)
    neg = np.array(neg)
    ax.scatter(*neg.T, c="red", zorder=10)
    plt.show()


def test_anapotcbm():
    calc_getter = AnaPotCBM
    # geom = AnaPotCBM().get_geom((0.818, 0.2233, 0.0))
    geom = AnaPotCBM().get_geom((0.2, 0.2, 0.0))
    # geom = AnaPotCBM().get_geom((0.3, 0.3, 0.0))
    # geom = AnaPotCBM().get_geom((0.4, 0.4, 0.0))
    # geom = AnaPotCBM().get_geom((0.5, 0.2, 0.0))
    geom = AnaPotCBM().get_geom((0.9, 0.8, 0.0))
    # geom = AnaPotCBM().get_geom((0.8, 0.7, 0.0))
    # geom = AnaPotCBM().get_geom((0.65, 0.7, 0.0))
    w, v = np.linalg.eigh(geom.hessian)
    print("eigenvals", w)
    N_imag = v[:,0]
    geoms = [geom, ]
    dimer_kwargs = {
        "ana_2dpot": True,
        "N_init": N_imag,
        "rot_type": "fourier",
        "trans_opt": "lbfgs",
        "f_tran_mod": False,
    }
    true_ts = (0, 0)
    dimer_result = dimer_method(geoms, calc_getter, **dimer_kwargs)
    dimer_cycles = dimer_result.dimer_cycles
    plot_dimer_cycles(dimer_cycles[-10:], pot=AnaPotCBM(), true_ts=true_ts)


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
        # "rot_f_thresh": 1e-4,
        # "rot_type": "direct",
        "f_thresh": 1e-4,
        "f_tran_mod": True,
        "multiple_translations": True,
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
        'mb_4': 16,
        'mb_5': 16,
    }
    results = dict()
    for to, tm in it.product(trans_opts, trans_memories):
        dimer_result = test_hcn_iso_dimer(to, tm)
        key = f"{to}_{tm}"
        results[key] = dimer_result.force_evals
        assert ref_evals[key] == dimer_result.force_evals
    print(results)


def test_baker_16():
    calc_kwargs = {
        "route": "HF/3-21G",
        "pal": 4,
        "mem": 1000,
        "charge": -1,
        "mult": 1,
    }
    def calc_getter():
        return Gaussian16(**calc_kwargs)
        # return XTB(**calc_kwargs)

    geom = geom_from_library("baker_ts/16_h2po4_anion.xyz")
    geom.set_calculator(calc_getter())
    geoms = [geom, ]

    downhill = geom_from_library("baker_ts/16_h2po4_anion_downhill.xyz")
    N_init = geom.coords - downhill.coords
    N_init /= np.linalg.norm(N_init)

    dimer_kwargs = {
        "max_step": 0.5,
        # 1e-2 Angstroem in bohr
        "dR_base": 0.0189,
        "rot_opt": "lbfgs",
        "trans_opt": "lbfgs",
        # "trans_opt": "mb",
        # "trans_memory": 10,
        "angle_tol": 5,
        "f_thresh": 1e-3,
        "max_cycles": 50,
        "f_tran_mod": True,
        "rot_type": "fourier",
        "multiple_translations": True,
        "N_init": N_init,
        "max_cycles": 15,
        "rot_f_thresh": 1.5e-3,
    }

    # 55 cycles; no multi trans, no f_tran_mod

    # import pdb; pdb.set_trace()
    dimer_result = dimer_method(geoms, calc_getter, **dimer_kwargs)
    return dimer_result

if __name__ == "__main__":
    # test_anapot_rot()
    # anapot_tester()
    # test_anapot_cbm_rot()
    # hcn_tester()
    plot_anapotcbm_curvature()
    # test_anapotcbm()
    # test_baker_16()
