#!/usr/bin/env python3

import numpy as np

from pysisyphus.intcoords import autodiff
from pysisyphus.helpers import geom_from_library


np.set_printoptions(suppress=True, precision=4)


def test_autodiff():
    geom = geom_from_library("h2o.xyz", coord_type="redund")
    # geom.jmol()
    print(geom)
    B_ref = geom.internal.B
    print(B_ref)

    sg = autodiff.stretch_grad
    bg = autodiff.bend_grad

    auto_funcs = {
        2: sg,
        3: bg,
    }

    int_ = geom.internal
    ref_funcs = {
        2: int_.calc_stretch,
        3: int_.calc_bend,
    }

    c3d = geom.coords3d
    for i, pc in enumerate(geom.internal._prim_internals):
        inds = pc.inds
        print(i, inds)
        l = len(inds)
        ag = auto_funcs[l](c3d, inds)
        _, rg = ref_funcs[l](c3d, inds, grad=True)
        np.testing.assert_allclose(ag.flatten(), rg)


if __name__ == "__main__":
    test_autodiff()
