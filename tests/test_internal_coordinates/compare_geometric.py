#!/usr/bin/env python3

import itertools as it

import numpy as np

from pysisyphus.helpers import geom_from_library


def h2o_comp():
    # geom = geom_from_library("h2o.xyz", coord_type="redund")
    geom = geom_from_library("h2o2_hf_321g_opt.xyz", coord_type="redund")
    # geom = geom_from_library("01_opt_aus_04_b3lyp_scan_min.xyz", coord_type="redund")
    xyz = geom.cart_coords
    int_ = geom.internal

    # from geometric.internal import Distance, Angle, Dihedral
    from internal import Distance, Angle, Dihedral
    from pysisyphus.intcoords.derivatives import d2q_b, d2q_a, d2q_d
    funcs = {
        2: Distance,
        3: Angle,
        4: Dihedral,
    }
    d2q_funcs = {
        2: d2q_b,
        3: d2q_a,
        4: d2q_d,
    }

    inds = it.chain(*int_.prim_indices)
    np.set_printoptions(precision=3, linewidth=120)
    c3d = geom.coords3d
    size_ = c3d.size
    for prim in int_._prim_internals[:2]:
        inds = prim.inds
        int_geom = funcs[len(inds)](*inds)
        val = int_geom.value(xyz)
        der = int_geom.derivative(xyz)
        print(int_geom, val)
        # if inds.size == 4:
            # import pdb; pdb.set_trace()
        np.testing.assert_allclose(der.flatten(), prim.grad)
        prim_coords = c3d[inds].flatten()
        d2q_ = d2q_funcs[len(inds)](*prim_coords)
        d2q_geom = int_geom.second_derivative(xyz)
        d2q = np.zeros_like(d2q_geom)
        cart_inds = list(it.chain(*[range(3*i,3*i+3) for i in inds]))
        flat_inds = [row*size_ + col for row, col in it.product(cart_inds, cart_inds)]
        print(cart_inds)
        print(flat_inds)
        print(size_)
        # import pdb; pdb.set_trace()
        # np.testing.assert_allclose(d2q, d2q_geom)
        # import pdb; pdb.set_trace()
        pass


if __name__ == "__main__":
    h2o_comp()
