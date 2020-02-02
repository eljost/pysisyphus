#!/usr/bin/env python3

# [1] https://doi.org/10.1016/S0022-2860(84)87198-7
#     Pulay, 1984
# [2] https://pubs.rsc.org/en/content/articlehtml/2002/cp/b108658h
#     Farkas, Schlegel, 2002


import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.XTB import XTB
from pysisyphus.helpers import geom_from_library
from pysisyphus.optimizers.gdiis import gdiis
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


def get_step(H, forces, trust):
    step = np.linalg.inv(H) @ forces
    norm = np.linalg.norm(step)
    if norm > trust:
        step = trust * step / norm
    return step


def test_sd_gdiis():
    geom = AnaPot.get_geom((0.4333, 3.14286, 0.))

    trust = 0.3
    H = np.eye(geom.coords.size)

    cs = list()
    dcs = list()
    fs = list()
    pairs = list()
    for i in range(50):
        forces = geom.forces

        cs.append(geom.coords)
        fs.append(forces)

        forces_norm = np.linalg.norm(forces)
        print(f"{i:02d}: {forces_norm:.6f}")
        if forces_norm < 1e-3:
            print("Converged!")
            break

        # Calculate reference step
        step = get_step(H, forces, trust)

        if len(fs) > 1:
            gdiis_kwargs = {
                "coords": cs,
                "forces": fs,
                "ref_step": step,
            }
            res = gdiis(fs, **gdiis_kwargs)
            if res:
                # Inter-/extrapolate coordinates and forces
                coords_ = res.coords
                forces = res.forces
                dcs.append(coords_)
                pairs.append((geom.coords.copy(), coords_, geom.coords + step))
                geom.coords = coords_
                forces_norm = np.linalg.norm(forces)
                # Get new step from DIIS coordinates & forces
                step = get_step(H, forces, trust)

        new_coords = geom.coords + step
        geom.coords = new_coords
    # return

    cs = np.array(cs)
    dcs = np.array(dcs)
    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    ax.plot(*cs.T[:2], "o-")
    for i, cs_ in enumerate(cs):
        ax.annotate(f"{i:02d}", cs_[:2])

    ax.plot(*dcs.T[:2], "o-")
    for oc, dc, rs in pairs:
        line = mlines.Line2D((oc[0], dc[0]), (oc[1], dc[1]), ls="--", color="k")
        ax.add_artist(line)
        line = mlines.Line2D((oc[0], rs[0]), (oc[1], rs[1]), ls="--", color="r")
        ax.add_artist(line)
    plt.show()


def test_rfo_benzene():
    geom = geom_from_library("benzene_shaken.xyz", coord_type="redund")
    calc = XTB(pal=4)
    geom.set_calculator(calc)
    opt = RFOptimizer(geom, dump=True)
    opt.run()


from pysisyphus.optimizers.gdiis import gediis

def from_coeffs(vec, coeffs):
    return np.sum(coeffs[:,None] * vec[::-1][:len(coeffs)], axis=0)


def test_sd_gediis():
    geom = AnaPot.get_geom((0.4333, 3.14286, 0.))

    trust = 0.3
    H = np.eye(geom.coords.size)

    ens = list()
    cs = list()
    dcs = list()
    fs = list()
    pairs = list()
    for i in range(25):
        forces = geom.forces
        energy = geom.energy

        cs.append(geom.coords)
        fs.append(forces)
        ens.append(energy)

        forces_norm = np.linalg.norm(forces)
        print(f"{i:02d}: {forces_norm:.6f}")
        if forces_norm < 1e-3:
            print("Converged!")
            break

        # Calculate reference step
        step = get_step(H, forces, trust)

        if len(fs) > 1:
            gediis_res = gediis(cs, ens, fs)


            if gediis_res is not None:
                # Inter-/extrapolate coordinates and forces
                coords_ = gediis_res.coords
                forces = gediis_res.forces
                dcs.append(coords_)
                geom.coords = coords_
                forces_norm = np.linalg.norm(forces)
                # Get new step from GEDIIS coordinates & forces
                step = get_step(H, forces, trust)

        new_coords = geom.coords + step
        geom.coords = new_coords

    cs = np.array(cs)
    dcs = np.array(dcs)
    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    ax.plot(*cs.T[:2], "o-")
    for i, cs_ in enumerate(cs):
        ax.annotate(f"{i:02d}", cs_[:2])

    ax.plot(*dcs.T[:2], "o-")
    for oc, dc, rs in pairs:
        line = mlines.Line2D((oc[0], dc[0]), (oc[1], dc[1]), ls="--", color="k")
        ax.add_artist(line)
        line = mlines.Line2D((oc[0], rs[0]), (oc[1], rs[1]), ls="--", color="r")
        ax.add_artist(line)
    plt.show()


def test_diis():
    geom = AnaPot.get_geom((0.4333, 3.14286, 0.))
    opt_kwargs = {
        # "max_step": 0.5,
        # "trust_max": 0.1,
        "thresh": "gau_tight",
        # "eins": True,
        "zwei": True,
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()
    calc = geom.calculator
    calc.plot()
    ax = calc.ax

    cs = np.array(opt.coords)
    ax.plot(*cs.T[:2], "ro-")
    plt.show()


def test_artemisin():
    # geom = geom_from_library("birkholz/artemisin.xyz", coord_type="redund")
    geom = geom_from_library("birkholz/artemisin.xyz")
    calc = XTB(charge=0, mult=1, pal=4)
    # geom = geom_from_library("birkholz/zn_edta.xyz", coord_type="redund")
    # geom = geom_from_library("birkholz/zn_edta.xyz")
    # calc = XTB(charge=-2, mult=1, pal=4)
    geom.set_calculator(calc)

    opt_kwargs_base = {
        "max_cycles": 50,
        # "max_cycles": 15,
        # "max_cycles": 8,
        "thresh": "gau",
        "trust_radius": 0.5,
        "trust_update": True,
        "hessian_update": "damped_bfgs",
        # "hessian_init": "fischer",
        "hessian_init": "calc",
        "hessian_recalc": 1,
        "line_search": True,
        # "hybrid": True,
        # "dump": True,
    }
    opt = RFOptimizer(geom, **opt_kwargs_base)
    opt.run()
    H = opt.H
    w, v = np.linalg.eigh(H)
    print(w)
    # import pdb; pdb.set_trace()


def test_hess_proj():
    geom = geom_from_library("birkholz/zn_edta.xyz")
    H = np.eye(geom.cart_coords.size)
    mwH = geom.mass_weigh_hessian(H)
    eH = geom.eckart_projection(mwH)
    w, v = np.linalg.eigh(eH)
    import pdb; pdb.set_trace()


def test_quad():
    def func(x):
        return x**2
    def dfunc(x):
        return 2*x
    coords = np.array(((-2, ), (1, )))
    energies = func(coords)
    forces = -dfunc(coords)
    gediis(coords, energies, forces)

if __name__ == "__main__":
    # test_sd_gdiis()
    # test_rfo_benzene()
    # test_sd_gediis()
    # test_diis()
    # test_artemisin()
    # test_hess_proj()
    test_quad()
