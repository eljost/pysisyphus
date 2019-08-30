#!/usr/bin/env python3

import os
from pathlib import Path

import numpy as np

from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.helpers import geom_from_xyz_file
from pysisyphus.run import do_final_hessian

THIS_DIR = Path(os.path.abspath(os.path.dirname(__file__)))

def test_do_final_hessian():
    fn = THIS_DIR  / "final_geometry.xyz"
    geom = geom_from_xyz_file(fn, coord_type="redund")
    calc = ORCA("")

    grad = np.load("grad.npy")
    hess = np.load("hess.npy")
    geom.hessian = hess
    geom.gradient = grad

    res = do_final_hessian(geom, save_hessian=False)

    neg_eigvals = res.neg_eigvals
    assert len(neg_eigvals) == 1
    np.testing.assert_allclose(neg_eigvals[0], -0.00224392407)


if __name__ == "__main__":
    test_do_final_hessian()
