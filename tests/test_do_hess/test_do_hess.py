#!/usr/bin/env python3

import numpy as np

from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.helpers import geom_from_xyz_file, do_final_hessian


def test_do_final_hessian(datadir):
    fn = datadir / "final_geometry.xyz"
    geom = geom_from_xyz_file(fn, coord_type="redund")
    calc = ORCA("")

    grad = np.load(datadir / "grad.npy")
    hess = np.load(datadir / "hess.npy")
    geom.hessian = hess
    geom.gradient = grad

    res = do_final_hessian(geom, save_hessian=False)

    neg_eigvals = res.neg_eigvals
    assert len(neg_eigvals) == 1
    np.testing.assert_allclose(neg_eigvals[0], -0.00224392407)
