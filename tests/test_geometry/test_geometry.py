#!/usr/bin/env python3

import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_xyz_file, geoms_from_trj
from pysisyphus.xyzloader import make_trj_str_from_geoms


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


def test_inertia_tensor():
    geom = geom_from_xyz_file("benzene.xyz")
    I = geom.inertia_tensor
    ref = (317.59537834, 317.59537834, 635.19075413)
    np.testing.assert_allclose(np.diag(I), ref)


def test_standard_orientation():
    geoms = geoms_from_trj("final.trj")[:50]
    aligned = [geom.std_orient for geom in geoms]
    with open("stdorient.trj", "w") as handle:
        handle.write(make_trj_str_from_geoms(geoms))


if __name__ == "__main__":
    # test_center_of_mass()
    # test_centroid()
    # test_inertia_tensor()
    test_standard_orientation()
