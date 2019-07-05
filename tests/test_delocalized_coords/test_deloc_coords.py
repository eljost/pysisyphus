#!/usr/bin/env python3

import numpy as np


from pysisyphus.calculators.XTB import XTB
from pysisyphus.calculators.Psi4 import Psi4
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.helpers import geom_from_library


np.set_printoptions(suppress=True, precision=6)


def test_dlc_constraints():
    xyz_fn = "fluorethylene.xyz"
    geom = geom_from_library(xyz_fn, coord_type="dlc")
    pi = (0, 1)
    i1 = geom.internal.get_index_of_prim_coord(pi)
    i2 = geom.internal.get_index_of_prim_coord((2,0,3))
    c1 = geom.internal.project_primitive_on_active_set(i1)
    c2 = geom.internal.project_primitive_on_active_set(i2)
    U = geom.internal.U
    V = np.concatenate((c1[:,None], c2[:,None], U), axis=1)
    b0 = V[0,:].copy()
    from pysisyphus.linalg import gram_schmidt
    gs = np.array(gram_schmidt(V))
    b0_ = gs[0,:].copy()
    import pdb; pdb.set_trace()


def run():
    # geom = geom_from_library("split.image_021.xyz")
    # geom.set_calculator(XTB())
    # opt = RFOptimizer(geom)
    # opt.run()

    # geom = geom_from_library("split.image_021.xyz", coord_type="redund")
    # geom.set_calculator(XTB())
    # opt = RFOptimizer(geom)
    # opt.run()

    xyz_fn = "fluorethylene.xyz"
    # xyz_fn = "split.image_021.xyz"
    # xyz_fn = "h2o2_hf_321g_opt.xyz"
    geom = geom_from_library(xyz_fn, coord_type="dlc")
    # geom = geom_from_library(xyz_fn, coord_type="redund")
    # int_ = geom.internal
    # int_.B
    # int_.set_delocalized_vectors()
    calc = XTB()
    # psi4_kwargs = {
        # "method": "hf",
        # "basis": "def2-svp",
    # }
    # calc = Psi4(**psi4_kwargs)
    geom.set_calculator(calc)
    opt_kwargs = {
        # "max_cycles": 1,
        "thresh": "gau_tight",
        "trust_max": 0.3,
        "trust_radius": 0.1,
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()


if __name__ == "__main__":
    # run()
    test_dlc_constraints()
