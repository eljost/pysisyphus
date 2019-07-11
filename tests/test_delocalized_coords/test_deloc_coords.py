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
    geom.internal.freeze_primitives(freeze)
    U = geom.internal.U
    assert U[0].sum() == approx(0)
    assert U[8].sum() == approx(0)


def run():
    xyz_fn = "fluorethylene.xyz"
    geom = geom_from_library(xyz_fn, coord_type="dlc")
    # freeze = ((0, 1), (2, 0, 3))
    # geom.internal.freeze_primitives(freeze)
    calc = XTB()
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
