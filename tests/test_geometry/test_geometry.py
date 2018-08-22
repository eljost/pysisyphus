#!/usr/bin/env python3

import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_xyz_file


def test_center_of_mass():
    geom = geom_from_xyz_file("benzene.xyz")
    R = geom.center_of_mass
    np.testing.assert_allclose(R, np.zeros(3), atol=1e-10)


if __name__ == "__main__":
    test_center_of_mass()
