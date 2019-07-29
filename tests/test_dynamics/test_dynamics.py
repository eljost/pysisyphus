#!/usr/bin/env python3

from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.MullerBrownSympyPot import MullerBrownPot
from pysisyphus.dynamics.velocity_verlet import md
from pysisyphus.dynamics.mdp import mdp


def test_velocity_verlet():
    geom = AnaPot.get_geom((0.52, 1.80, 0))
    x0 = geom.coords.copy()
    v0 = .1 * np.random.rand(*geom.coords.shape)
    t = 25
    dts = (.005, .01, .02, .04, .08)
    all_xs = list()
    for dt in dts:
        geom.coords = x0.copy()
        md_kwargs = {
            "v0": v0.copy(),
            "t": t,
            "dt": dt,
        }
        md_result = md(geom, **md_kwargs)
        all_xs.append(md_result.coords)
    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    for dt, xs in zip(dts, all_xs):
        # ax.plot(*xs.T[:2], "o-", label=f"dt={dt:.3f}")
        ax.plot(*xs.T[:2], "-", label=f"dt={dt:.3f}")
    ax.legend()
    plt.show()


def test_mdp():
    coords = (-0.82200156,  0.6243128, 0)
    geom = MullerBrownPot.get_geom(coords)

    A = (-0.5592, 1.443, 0)
    B = (0.605, 0.036, 0)
    rad = 0.05
    def stopA(x, rad=rad):
        return np.linalg.norm(x-A) < rad
    def stopB(x, rad=rad):
        return np.linalg.norm(x-B) < rad
    term_funcs = (stopA, stopB)

    mdp_kwargs = {
        "E_excess": 0.1,
        "term_funcs": term_funcs,
        "epsilon": 5e-4,
        "ascent_alpha": 0.05,
        "t_init": 0.15,
        "t": 3,
        "dt": 0.001,
    }
    # np.random.seed(25032018)
    res = mdp(geom, **mdp_kwargs)

    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    ax.plot(*res.ascent_xs.T[:2], "ro-")
    ax.plot(*res.md_init_plus.coords.T[:2], "-", lw=3)
    ax.plot(*res.md_init_minus.coords.T[:2], "-", lw=3)
    cA = Circle(A[:2], radius=rad)
    ax.add_artist(cA)
    cB = Circle(B[:2], radius=rad)
    ax.add_artist(cB)
    ax.plot(*res.md_fin_plus.coords.T[:2], "-", lw=3)
    ax.plot(*res.md_fin_minus.coords.T[:2], "-", lw=3)

    plt.show()




if __name__ == "__main__":
    test_velocity_verlet()
    test_mdp()
