from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.MullerBrownSympyPot import MullerBrownPot
from pysisyphus.dynamics import mdp


def test_muller_brown_mdp():
    coords = (-0.82200156,  0.6243128, 0)
    geom = MullerBrownPot.get_geom(coords)

    # Termination functions
    A = (-0.5592, 1.443, 0)
    B = (0.605, 0.036, 0)
    radius = 0.05
    def stopA(x, rad=radius):
        return np.linalg.norm(x-A) < radius
    def stopB(x, rad=radius):
        return np.linalg.norm(x-B) < radius
    term_funcs = {
        "nearA": stopA,
        "nearB": stopB,
    }

    mdp_kwargs = {
        "E_excess": 0.2,
        "term_funcs": term_funcs,
        "epsilon": 5e-4,
        "ascent_alpha": 0.0125,
        "t_init": 0.15,
        "t": 3,
        "dt": 0.001,
    }
    np.random.seed(25032018)
    res = mdp(geom, **mdp_kwargs)

    assert res.md_fin_plus_term == "nearB"
    assert res.md_fin_minus_term == "nearA"

    # calc = geom.calculator
    # calc.plot()
    # ax = calc.ax
    # ax.plot(*res.ascent_xs.T[:2], "ro-")
    # ax.plot(*res.md_init_plus.coords.T[:2], "-", lw=3)
    # ax.plot(*res.md_init_minus.coords.T[:2], "-", lw=3)
    # cA = Circle(A[:2], radius=radius)
    # ax.add_artist(cA)
    # cB = Circle(B[:2], radius=radius)
    # ax.add_artist(cB)
    # ax.plot(*res.md_fin_plus.coords.T[:2], "-", lw=3)
    # ax.plot(*res.md_fin_minus.coords.T[:2], "-", lw=3)
    # plt.show()
