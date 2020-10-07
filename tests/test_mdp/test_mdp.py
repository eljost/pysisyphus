from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import numpy as np
import pytest

from pysisyphus.calculators.MullerBrownSympyPot import MullerBrownPot
from pysisyphus.calculators import XTB
from pysisyphus.dynamics import mdp
from pysisyphus.dynamics.mdp import parse_raw_term_funcs
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using


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

@using("xtb")
@pytest.mark.parametrize(
    "external_md", [
        False,
        # True,  # term_funcs not implemented
    ]
)
def test_so3hcl_diss_mdp(external_md):
    """See
        https://aip.scitation.org/doi/pdf/10.1063/1.5082885
    """
    geom = geom_loader("lib:so3hcl_diss_ts_opt.xyz")
    geom.set_calculator(XTB(pal=4))
    # geom.jmol()

    def hcl_dissociated(coords):
        coords3d = coords.reshape(-1, 3)
        chlorine = coords3d[5]
        sulphur = coords3d[0]
        return np.linalg.norm(chlorine-sulphur) > 4.6  # Bohr (> 2.4 Å)

    def oh_formed(coords):
        coords3d = coords.reshape(-1, 3)
        oxygen = coords3d[1]
        hydrogen = coords3d[4]
        return np.linalg.norm(oxygen-hydrogen) < 1.9  # Bohr (< Å)

    term_funcs = {
        "hcl_dissociated": hcl_dissociated,
        "oh_formed": oh_formed,
    }

    mdp_kwargs = {
        # About 5 kcal/mol as given in the paper
        "E_excess": 0.0079,
        "epsilon": 5e-4,
        "ascent_alpha": 0.025,
        "term_funcs": term_funcs,
        "t_init": 20,
        # Paper uses 200
        "t": 100,
        "dt": .5,
        "seed": 9081302,
        "external_md": external_md,
        "max_init_trajs": 1,
    }
    res = mdp(geom, **mdp_kwargs)

    assert res.md_fin_plus_term == "hcl_dissociated"
    assert res.md_fin_minus_term == "oh_formed"


@pytest.mark.parametrize(
    "distance, lt_ref, gt_ref, le_ref, ge_ref, eq_ref",
    [
        (1.0,  True, False,  True, False, False),
        (2.0,  True, False,  True, False, False),
        (2.9,  True, False,  True, False, False),
        (3.0, False, False,  True,  True, True),
        (3.1, False,  True, False,  True, False),
        (3.2, False,  True, False,  True, False),
    ],
)
def test_term_func_dissociation(distance, lt_ref, gt_ref, le_ref, ge_ref, eq_ref):
    raw_term_funcs = {
        "lt": "0,1 < 3.0",
        "gt": "0,1 > 3.0",
        "le": "0,1 <= 3.0",
        "ge": "0,1 >= 3.0",
        "eq": "0,1 == 3.0",
    }
    term_funcs = parse_raw_term_funcs(raw_term_funcs)
    coords3d = np.zeros((2, 3))
    coords3d[1,2] = distance
    keys = ("lt", "gt", "le", "ge", "eq")
    lt_res, gt_res, le_res, ge_res, eq_res = [term_funcs[key](coords3d) for key in keys]

    assert lt_res == lt_ref
    assert gt_res == gt_ref
    assert le_res == le_ref
    assert ge_res == ge_ref
    assert eq_res == eq_ref
