#!/usr/bin/env python3

import itertools as it
from pathlib import Path

import numpy as np

from pysisyphus.helpers import geom_from_library
from pysisyphus.calculators.Gaussian16 import Gaussian16

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
    geom = geom_from_library(xyz_fn, coord_type="redund")

    calc_kwargs = {
        "route": "HF/3-21G",
        "charge": 0,
        "mult": 1,
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


if __name__ == "__main__":
    bond_ref()
    test_h2()
    # test_h2o()
    # test_h2o2()
