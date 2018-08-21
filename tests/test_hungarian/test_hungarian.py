#!/usr/bin/env python3

import numpy as np

from pysisyphus.helpers import geom_from_xyz_file, match_geoms


def get_geoms():
    benz = geom_from_xyz_file("benzene.xyz")
    # Here the atom at index 3 was moved to index 5
    benz_mod = geom_from_xyz_file("benzene_mod.xyz")

    return benz, benz_mod


def test_benzene_swap():
    benz, benz_mod = get_geoms()
    match_geoms(benz, benz_mod)
    np.testing.assert_allclose(benz.coords3d, benz_mod.coords3d)


def test_benzene_swap2():
    benz, benz_mod = get_geoms()
    match_geoms(benz_mod, benz)
    np.testing.assert_allclose(benz.coords3d, benz_mod.coords3d)


if __name__ == "__main__":
    test_benzene_swap()
