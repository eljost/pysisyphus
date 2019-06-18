#!/usr/bin/env python3

import itertools as it

import numpy as np

from pysisyphus.helpers import geom_from_library


def h2o_comp():
    # geom = geom_from_library("h2o.xyz", coord_type="redund")
    # geom = geom_from_library("h2o2_hf_321g_opt.xyz", coord_type="redund")
    geom = geom_from_library("01_opt_aus_04_b3lyp_scan_min.xyz", coord_type="redund")
    xyz = geom.cart_coords
    int_ = geom.internal

    from geometric.internal import Distance, Angle, Dihedral
    funcs = {
        2: Distance,
        3: Angle,
        4: Dihedral,
    }

    inds = it.chain(*int_.prim_indices)
    np.set_printoptions(precision=3, linewidth=120)
    for prim in int_._prim_internals:
        inds = prim.inds
        int_geom = funcs[len(inds)](*inds)
        val = int_geom.value(xyz)
        der = int_geom.derivative(xyz)
        print(int_geom, val)
        np.testing.assert_allclose(der.flatten(), prim.grad)


if __name__ == "__main__":
    h2o_comp()
