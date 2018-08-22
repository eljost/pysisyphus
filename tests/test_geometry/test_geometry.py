#!/usr/bin/env python3

import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_xyz_file


def test_center_of_mass():
    geom = geom_from_xyz_file("benzene.xyz")
    R = geom.center_of_mass
    np.testing.assert_allclose(R, np.zeros(3), atol=1e-10)


def test_centroid():
    geom = geom_from_xyz_file("benzene.xyz")
    ref = np.zeros(3)
    np.testing.assert_allclose(geom.centroid, ref, atol=1e-10)
    c3d = geom.coords3d
    x_translate = (1, 0, 0)
    c3d += x_translate
    ref += x_translate
    np.testing.assert_allclose(geom.centroid, ref, atol=1e-10)


if __name__ == "__main__":
    test_center_of_mass()
    test_centroid()
