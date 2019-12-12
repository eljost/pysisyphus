#!/usr/bin/env python3

# [1] https://aip.scitation.org/doi/pdf/10.1063/1.1515483

import itertools as it

import numpy as np

from pysisyphus.helpers import geom_from_library
from pysisyphus.calculators.XTB import XTB
from pysisyphus.intcoords.derivatives import dq_b, d2q_b, dq_a, d2q_a, dq_d, d2q_d


def test_wilson():
    geom = geom_from_library("split.image_021.xyz", coord_type="redund")
    forces_fn = "forces"
    hessian_fn = "hessian"
    # geom.set_calculator(XTB(pal=4))
    # forces = geom.forces
    # _ = geom.hessian
    # hessian = geom._hessian
    # np.savetxt(forces_fn, forces)
    # np.savetxt(hessian_fn, hessian)
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
    for i, pc in enumerate(int_._prim_internals):
        inds = pc.inds
        print(f"{i:02d}: {v[len(inds)]}")
        g = pc.grad
        g_ = grad_wrapper(inds)
        try:
            np.testing.assert_allclose(g, g_)
        except:
            np.testing.assert_allclose(g, -g_)

    dg_funcs = {
        2: d2q_b,
        3: d2q_a,
        4: d2q_d,
    }
    def grad_deriv_wrapper(inds):
        coords_flat = c3d[inds].flatten()
        dgrad = dg_funcs[len(inds)](*coords_flat)
        return dgrad

    gradient = -np.loadtxt(forces_fn)
    hessian = np.loadtxt(hessian_fn)
    size_ = geom.cart_coords.size
    K_flat = np.zeros(size_ * size_)

    assert len(gradient) == len(int_._prim_internals)
    for pc, g in zip(int_._prim_internals, gradient):
        # Contract with gradient
        dg = g * grad_deriv_wrapper(pc.inds)
        # Depending on the type of internal coordinate dg is a flat array
        # of size 36 (stretch), 81 (bend) or 144 (torsion).
        #
        # An internal coordinate contributes to an element K[j, k] of the
        # K matrix if the cartesian coordinate indices j and k belong to an
        # atom that contributes to the respective internal coordinate.
        #
        # As for now we build up the K matrix as flat array. To add the dg
        # entries at the appropriate places in K_flat we have to calculate
        # the corresponding flat indices of dg in K_flat.
        cart_inds = list(it.chain(*[range(3*i,3*i+3) for i in pc.inds]))
        flat_inds = [row*size_ + col for row, col in it.product(cart_inds, cart_inds)]
        K_flat[flat_inds] += dg
    K = K_flat.reshape(-1, size_)

    # K[np.abs(K) < 1e-10] = np.nan
    # np.testing.assert_allclose(K,K.T)

    # Transform hessian
    H_int_ = int_.transform_hessian(hessian)
    H_int = int_.transform_hessian(hessian - K)
    import pdb; pdb.set_trace()
    pass


if __name__ == "__main__":
    test_wilson()
