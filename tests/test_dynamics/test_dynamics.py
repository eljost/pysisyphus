#!/usr/bin/env python3

# [1] https://doi.org/10.1063/1.5082885

from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.constants import BOHR2ANG
from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.MullerBrownSympyPot import MullerBrownPot
from pysisyphus.calculators.XTB import XTB
from pysisyphus.dynamics.helpers import get_velocities
from pysisyphus.dynamics.mdp import mdp
from pysisyphus.dynamics.velocity_verlet import md
from pysisyphus.helpers import geom_from_library


def test_velocity_verlet():
    geom = AnaPot.get_geom((0.52, 1.80, 0))
    x0 = geom.coords.copy()
    v0 = .1 * np.random.rand(*geom.coords.shape)
    t = 3
    dts = (.005, .01, .02, .04, .08)
    all_xs = list()
    for dt in dts:
        geom.coords = x0.copy()
        md_kwargs = {
            "v0": v0.copy(),
            "t": t,
            "dt": dt,
            "vcom": True,
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


def test_so3hcl_diss():
    """See [1]"""
    def get_geom():
        geom = geom_from_library("so3hcl_diss_ts_opt.xyz")
        geom.set_calculator(XTB(pal=4))
        return geom

    geom = get_geom()
    mdp_kwargs = {
        # About 5 kcal/mol
        "E_excess": 0.0079,
        "term_funcs": list(),
        "epsilon": 5e-4,
        "ascent_alpha": 0.05,
        "t_init": 20,
        # Paper uses 200
        "t": 100,
        "dt": .5,
        "seed": 25032018,
        # "external_md": True,
        "max_init_trajs": 1,
    }
    res = mdp(geom, **mdp_kwargs)

    # geom = get_geom()
    # mdp_kwargs["E_excess"] = 0
    # res_ee = mdp(geom, **mdp_kwargs)


def test_so3hcl_md():
    geom = geom_from_library("so3hcl_diss_ts_opt.xyz")
    geom.set_calculator(XTB(pal=4))

    v0 = .025 * np.random.rand(*geom.coords.shape)
    md_kwargs = {
        "v0": v0,
        "t": 400,
        "dt": 1,
    }
    res = md(geom, **md_kwargs)

    from pysisyphus.xyzloader import make_trj_str
    def dump_coords(coords, trj_fn):
        coords = np.array(coords)
        coords = coords.reshape(-1, len(geom.atoms), 3) * BOHR2ANG
        trj_str = make_trj_str(geom.atoms, coords)
        with open(trj_fn, "w") as handle:
            handle.write(trj_str)
    dump_coords(res.coords, "md.trj")


def test_xtb_md():
    geom = geom_from_library("so3hcl_diss_ts_opt.xyz")
    calc = XTB(pal=4)

    T = 298.15
    velocities = get_velocities(geom, T=T)
    geoms = calc.run_md(geom.atoms, geom.coords, t=200, step=0.1,
                        velocities=velocities)


def test_oniom_md():
    calc_dict = {
        "high": {
            "type": "pypsi4",
            "method": "scf",
            "basis": "sto-3g",
        },
        "low": {
            "type": "pyxtb",
        },
    }
    high_inds = (4,5,6)
    from pysisyphus.calculators.ONIOM import ONIOM
    oniom = ONIOM(calc_dict, high_inds)

    geom = geom_from_library("acetaldehyd_oniom.xyz")
    geom.set_calculator(oniom)

    v0 = .005 * np.random.rand(*geom.coords.shape)
    md_kwargs = {
        "v0": v0,
        "t": 40,
        "dt": 0.5,
    }
    md_result = md(geom, **md_kwargs)
    from pysisyphus.xyzloader import make_trj_str

    coords = md_result.coords.reshape(-1, len(geom.atoms), 3) * BOHR2ANG
    trj_str = make_trj_str(geom.atoms, coords)
    with open("md.trj", "w") as handle:
        handle.write(trj_str)
    # import pdb; pdb.set_trace()


if __name__ == "__main__":
    test_velocity_verlet()
    # test_mdp()
    # test_so3hcl_diss()
    # test_so3hcl_md()
    # test_xtb_md()
    # test_oniom_md()
