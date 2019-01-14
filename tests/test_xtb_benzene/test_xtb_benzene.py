#!/usr/bin/env python3

import os
from pathlib import Path

import pytest
import numpy as np

from pysisyphus.calculators.XTB import XTB
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_library

from qchelper.geometry import parse_xyz_file

THIS_DIR = Path(os.path.dirname(__file__))


@pytest.mark.skip
def test_gradient():
    geom = geom_from_library("benzene.xyz")
    geom.set_calculator(XTB(gfn=1))
    grad = geom.gradient


    reference_fn = THIS_DIR / "gradient.reference"
    reference = np.loadtxt(reference_fn)

    # Quite some difference between xTB 5.8 and 4.9.4 ...
    np.testing.assert_allclose(reference.flatten(), grad, rtol=1e-5)


if __name__ == "__main__":
    test_gradient()
