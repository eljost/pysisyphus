#!/usr/bin/env python3

import os
from pathlib import Path

import numpy as np
import rmsd

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_xyz_file, match_geoms
from pysisyphus.stocastic.align import match_geom_atoms

THIS_DIR = Path(os.path.dirname(os.path.realpath(__file__)))


def get_geoms():
    benz = geom_from_xyz_file(THIS_DIR / "benzene.xyz")
    # Here the atom at index 3 was moved to index 5
    benz_mod = geom_from_xyz_file(THIS_DIR / "benzene_mod.xyz")

    return benz, benz_mod


def test_benzene_swap():
    benz, benz_mod = get_geoms()
    match_geoms(benz, benz_mod)
    np.testing.assert_allclose(benz.coords3d, benz_mod.coords3d)


def test_benzene_swap2():
    benz, benz_mod = get_geoms()
    match_geoms(benz_mod, benz)
    np.testing.assert_allclose(benz.coords3d, benz_mod.coords3d)


def shuffle_geometry(geom):
    """Shuffles the coordinates of a certain atom type.

    Used in testing the match_geom_atoms function.

    Returns
    -------
    shuffled_geom : Geometry
        Geometry containing shuffled coordinates.
    """
    coords_dict, inds = geom.coords_by_type
    shuffled_geom = Geometry(geom.atoms, geom.coords.copy())
    c3d = shuffled_geom.coords3d
    for atom, coords in coords_dict.items():
        np.random.shuffle(coords)
        c3d[inds[atom]] = coords
    return shuffled_geom


def test_arbitrary_shuffles():
    """Given an original geometry and a shuffled geometry derived
    from the original one, the Hungarian method should be able to
    fully reconstruct the original order of atoms."""

    benz, _ = get_geoms()
    coords_by_type, inds_by_type = benz.coords_by_type
    atom_num = len(benz.atoms)
    initial_order = np.arange(atom_num)
    atom_arr = np.array(benz.atoms)
    c3d = benz.coords3d.copy()

    for i in range(10):
        shuffled_geom = shuffle_geometry(benz)

        org_rmsd = rmsd.rmsd(benz.coords3d, shuffled_geom.coords3d)
        matched_geom = match_geom_atoms(benz, shuffled_geom, hydrogen=True)
        matched_rmsd = rmsd.rmsd(benz.coords3d, matched_geom.coords3d)
        print(f"RMSD(shuffled) {org_rmsd:.2f} "
              f"RMSD(matched) {matched_rmsd:.2f}")
        assert np.allclose(matched_rmsd, 0.0)


if __name__ == "__main__":
    # test_benzene_swap()
    test_arbitrary_shuffles()
