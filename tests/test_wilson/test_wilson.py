#!/usr/bin/env python3

# [1] https://aip.scitation.org/doi/pdf/10.1063/1.1515483

import numpy as np

from pysisyphus.helpers import geom_from_library
from pysisyphus.calculators.XTB import XTB
from wilson_gen import dq_b, d2q_b, dq_a, d2q_a, dq_d, d2q_d


def test_wilson():
    geom = geom_from_library("split.image_021.xyz", coord_type="redund")
    fn = "forces"
    # geom.set_calculator(XTB(pal=4))
    # forces = geom.forces
    # np.savetxt(fn, forces)
    int_ = geom.internal
    c3d = geom.coords3d

    g_funcs = {
        2: dq_b,
        3: dq_a,
        4: dq_d,
    }
    zero_row = np.zeros_like(geom.coords3d)
    def grad_wrapper(inds):
        coords_flat = c3d[inds].flatten()
        grad = g_funcs[len(inds)](*coords_flat)
        B_row = zero_row.copy()
        B_row[inds] = grad.reshape(-1, 3)
        return B_row.flatten()

    v = {
        2: "stretch",
        3: "bend",
        4: "torsion",
    }
    for i, pc in enumerate(int_._prim_coords):
        inds = pc.inds
        print(f"{i:02d}: {v[len(inds)]}")
        g = pc.grad
        g_ = grad_wrapper(inds)
        try:
            np.testing.assert_allclose(g, g_)
        except:
            np.testing.assert_allclose(g, -g_)

    # def s(stre, bend, tors):
        # return 8*(stre*36 + bend*81 + tors*144) / 1e6

    dg_funcs = {
        2: d2q_b,
        3: d2q_a,
        4: d2q_d,
    }
    def grad_deriv_wrapper(inds):
        coords_flat = c3d[inds].flatten()
        dgrad = dg_funcs[len(inds)](*coords_flat)
        return dgrad

    # dgs = [grad_deriv_wrapper(pc.inds) for pc in int_._prim_coords]

    grad = -np.loadtxt(fn)
    # Contract values with internal gradient, Eq. (8) in [1]
    dgs = [g*grad_deriv_wrapper(pc.inds) for pc, g in zip(int_._prim_coords, grad)]
    import pdb; pdb.set_trace()
    pass
    # for i, pc in enumerate(int_._prim_coords):
        # inds = pc.inds
        # dg = grad_deriv_wrapper(pc.inds)
        # import pdb; pdb.set_trace()
        # pass


if __name__ == "__main__":
    test_wilson()
