#!/usr/bin/env python3

import itertools as it
from pathlib import Path

import numpy as np

from pysisyphus.helpers import geom_from_library, geom_from_xyz_file
from pysisyphus.calculators.Gaussian16 import Gaussian16
from pysisyphus.calculators.Psi4 import Psi4

np.set_printoptions(precision=2, linewidth=120)


def bond_ref():
    xyz_fn = "h2.xyz"
    geom = geom_from_library(xyz_fn, coord_type="redund")
    def d2(geom):
        c3d = geom.coords3d
        n, m = 0, 1
        N_ = c3d[n]
        M_ = c3d[m]
        U_ = M_ - N_
        u_norm = np.linalg.norm(U_)
        U = U_ / u_norm
        inds = (n, m)

        derivs = list()
        kr = lambda i, j: 1 if i == j else 0
        for a, i in it.product((0, 1), (0, 1, 2)):
            for b, j in it.product((0, 1), (0, 1, 2)):
                d = (-1)**kr(a, b) * (U[i]*U[j] - kr(i, j))/u_norm
                print(a, b, i,  j, d)
                derivs.append(d)
        return derivs
    # import pdb; pdb.set_trace()
    d = d2(geom)
    d = np.array(d)
    d = d.reshape(-1, 6)
    import pdb; pdb.set_trace()
    np.testing.assert_allclose(d, d.T)
    import pdb; pdb.set_trace()
    pass



def test_base(xyz_fn):
    try:
        geom = geom_from_library(xyz_fn, coord_type="redund")
    except FileNotFoundError:
        geom = geom_from_xyz_file(xyz_fn, coord_type="redund")

    import pdb; pdb.set_trace()
    calc_kwargs = {
        # "route": "HF/3-21G",
        "route": "B3LYP/6-31G*",
        "charge": 0,
        "mult": 1,
        "pal": 4,
    }
    calc = Gaussian16(**calc_kwargs)
    geom.set_calculator(calc)

    H = geom.hessian
    xyz_path = Path(xyz_fn)
    H_fn = xyz_path.with_suffix(".hessian")
    np.savetxt(H_fn, H)
    print(f"Wrote hessian of '{xyz_fn}' to '{H_fn}'")


def test_h2():
    xyz_fn = "h2.xyz"
    test_base(xyz_fn)


def test_h2o():
    xyz_fn = "h2o.xyz"
    import pdb; pdb.set_trace()
    test_base(xyz_fn)


def test_h2o2():
    xyz_fn = "h2o2_short.xyz"
    test_base(xyz_fn)


def test_bz():
    xyz_fn = "benzene_ref.xyz"
    test_base(xyz_fn)


def test_dopamine():
    xyz_fn = "dopamine.xyz"
    test_base(xyz_fn)


def test_h2o2_opt():
    xyz_fn = "h2o2_hf_def2svp_opt.xyz"
    geom = geom_from_xyz_file(xyz_fn, coord_type="redund")

    cart_H_fn = "h2o2_opt_cart.hessian"
    xyz_path = Path(xyz_fn)
    H_fn = xyz_path.with_suffix(".hessian")
    cart_forces_fn = "h2o2_opt_cart.forces"


    calc_kwargs = {
        "method": "HF",
        "basis": "def2-svp",
        "charge": 0,
        "mult": 1,
        "pal": 4,
    }
    calc = Psi4(**calc_kwargs)
    geom.set_calculator(calc)

    cart_H = np.loadtxt(cart_H_fn)
    cart_f = np.loadtxt(cart_forces_fn)
    geom._forces = cart_f
    geom._hessian = cart_H
    import pdb; pdb.set_trace()

    H = geom.hessian
    np.savetxt(H_fn, H)
    np.savetxt(cart_H_fn, geom._hessian)
    np.savetxt(cart_forces_fn, geom._forces)
    print(f"Wrote hessian of '{xyz_fn}' to '{H_fn}'")
    # import pdb; pdb.set_trace()

if __name__ == "__main__":
    # bond_ref()
    # test_h2()
    # test_h2o()
    # test_h2o2()
    # test_bz()
    # test_dopamine()
    test_h2o2_opt()
