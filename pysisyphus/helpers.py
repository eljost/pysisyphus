import os
from pathlib import Path

import numpy as np
import scipy as sp

from pysisyphus.constants import ANG2BOHR
from pysisyphus.Geometry import Geometry
from pysisyphus.xyzloader import parse_xyz_file, parse_trj_file


def geom_from_xyz_file(xyz_fn, **kwargs):
    atoms, coords = parse_xyz_file(xyz_fn)
    coords *= ANG2BOHR
    geom = Geometry(atoms, coords.flatten(), **kwargs)
    return geom


def geom_from_library(xyz_fn, **kwargs):
    this_dir = os.path.abspath(os.path.dirname(__file__))
    xyz_dir = Path(this_dir) / "../xyz_files/"
    atoms, coords = parse_xyz_file(xyz_dir / xyz_fn)
    coords *= ANG2BOHR
    geom = Geometry(atoms, coords.flatten(), **kwargs)
    return geom


def geoms_from_trj(trj_fn, **kwargs):
    geoms = [Geometry(atoms, coords.flatten()*ANG2BOHR, **kwargs)
             for atoms, coords in parse_trj_file(trj_fn)
    ]
    return geoms


def procrustes(geometry):
    # http://nghiaho.com/?page_id=671#comment-559906
    image0 = geometry.images[0]
    coords3d = image0.coords.reshape(-1, 3)
    centroid = coords3d.mean(axis=0)
    last_centered = coords3d - centroid
    geometry.set_coords_at(0, last_centered.flatten())
    atoms_per_image = len(image0.atoms)

    # Don't rotate the first image, so just add identitiy matrices
    # for every atom.
    rot_mats = [np.eye(3)]*atoms_per_image
    for i, image in enumerate(geometry.images[1:], 1):
        coords3d = image.coords.reshape(-1, 3)
        centroid = coords3d.mean(axis=0)
        # Center next image
        centered = coords3d - centroid
        tmp_mat = centered.T.dot(last_centered)
        U, W, Vt = np.linalg.svd(tmp_mat)
        rot_mat = U.dot(Vt)
        # Avoid reflections
        if np.linalg.det(rot_mat) < 0:
            U[:, -1] *= -1
            rot_mat = U.dot(Vt)
        # Rotate the coords
        rotated3d = centered.dot(rot_mat)
        geometry.set_coords_at(i, rotated3d.flatten())
        last_centered = rotated3d
        # Append one rotation matrix per atom
        rot_mats.extend([rot_mat]*atoms_per_image)
    return rot_mats


def fit_rigid(geometry, vectors=(), hessian=None):
    rot_mats = procrustes(geometry)
    G = sp.linalg.block_diag(*rot_mats)
    rotated_vectors = [vec.dot(G) for vec in vectors]
    if hessian is None:
        return rotated_vectors

    #rotated_hessian = G.dot(hessian).dot(G.T)
    rotated_hessian = G.T.dot(hessian).dot(G)
    #rotated_hessian = G*hessian*G.T
    return rotated_vectors, rotated_hessian


if __name__ == "__main__":
    print(load_geometry("hcn.xyz"))
    print(geoms_from_trj("cycle_040.trj"))
