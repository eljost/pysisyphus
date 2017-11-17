#!/usr/bin/env python3

import os
from pathlib import Path

import numpy as np

from pysisyphus.calculators.XTB import XTB
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_library

from qchelper.geometry import parse_xyz_file

THIS_DIR = Path(os.path.dirname(__file__))


def test_gradient():
    geom = geom_from_library("benzene.xyz")
    geom.set_calculator(XTB())
    grad = geom.gradient


    reference_fn = THIS_DIR / "gradient.reference"
    reference = np.loadtxt(reference_fn)

    np.testing.assert_allclose(reference.flatten(), grad)


if __name__ == "__main__":
    test_gradient()
