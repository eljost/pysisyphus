#!/usr/bin/env python3

import numpy as np
from pytest import approx


from pysisyphus.calculators.XTB import XTB
from pysisyphus.calculators.Psi4 import Psi4
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.helpers import geom_from_library


np.set_printoptions(suppress=True, precision=3, linewidth=150)


def test_dlc_constraints():
    xyz_fn = "fluorethylene.xyz"
    geom = geom_from_library(xyz_fn, coord_type="dlc")
    freeze = ((0, 1), (2, 0, 3))
    int_ = geom.internal
    int_.freeze_primitives(freeze)
    U = int_.U
    constraints = int_.constraints
    assert U[0,len(freeze):].sum() == approx(0)
    assert U[8,len(freeze):].sum() == approx(0)
    int_.reset_constraints()
    constraints_ = int_.constraints


def test_h2o2_dlc_constraints():
    xyz_fn = "h2o2_hf_321g_opt.xyz"
    geom = geom_from_library(xyz_fn, coord_type="dlc")
    xyz_fn2 = "h2o2_rot2.xyz"
    geom2 = geom_from_library(xyz_fn2, coord_type="dlc")
    tangent = geom2 - geom
    from pysisyphus.intcoords.helpers import get_tangent
    prim_tangent = get_tangent(geom2.internal.prim_coords,
                               geom.internal.prim_coords,
                               geom2.internal.dihed_start)

    freeze = ((0, 1, 2, 3), )
    int_ = geom.internal
    int_.freeze_primitives(freeze)
    U = int_.U

    constraints = int_.constraints
    Uc = np.concatenate((int_.constraints, int_.U), axis=1)
    Uc_inv = np.linalg.pinv(Uc.T)
    dlc_tangent = Uc.T.dot(prim_tangent/np.linalg.norm(prim_tangent))
    # int_.reset_constraints()
    # constraints_ = int_.constraints
    import pdb; pdb.set_trace()


def run():
    xyz_fn = "fluorethylene.xyz"
    geom = geom_from_library(xyz_fn, coord_type="dlc")
    freeze = ((0, 1), (2, 0, 3))
    geom.internal.freeze_primitives(freeze)

    # XTB
    calc = XTB()
    # Psi4
    # from pysisyphus.calculators.Psi4 import Psi4
    # calc = Psi4(method="hf", basis="sto-3g")

    geom.set_calculator(calc)
    opt_kwargs = {
        # "max_cycles": 1,
        "thresh": "gau_tight",
        "trust_max": 0.3,
        "trust_radius": 0.1,
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()
    assert opt.is_converged


if __name__ == "__main__":
    run()
    test_dlc_constraints()
    test_h2o2_dlc_constraints()
