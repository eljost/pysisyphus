import os
from pathlib import Path
import re

import numpy as np
import scipy as sp
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist

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


def geoms_from_trj(trj_fn, first=None, **kwargs):
    geoms = [Geometry(atoms, coords.flatten()*ANG2BOHR, **kwargs)
             for atoms, coords in parse_trj_file(trj_fn)[:first]
    ]
    return geoms


def align_geoms(geoms):
    # http://nghiaho.com/?page_id=671#comment-559906
    first_geom = geoms[0]
    coords3d = first_geom.coords.reshape(-1, 3)
    centroid = coords3d.mean(axis=0)
    last_centered = coords3d - centroid
    first_geom.coords3d = last_centered
    atoms_per_image = len(first_geom.atoms)

    # Don't rotate the first image, so just add identitiy matrices
    # for every atom.
    rot_mats = [np.eye(3)]*atoms_per_image
    for i, geom in enumerate(geoms[1:], 1):
        coords3d = geom.coords3d
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
        geom.coords3d = rotated3d
        last_centered = rotated3d
        # Append one rotation matrix per atom
        rot_mats.extend([rot_mat]*atoms_per_image)
    return rot_mats


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


def chunks(l, n):
    """Yield successive n-sized chunks from l.
    https://stackoverflow.com/a/312464
    """
    for i in range(0, len(l), n):
        yield l[i:i + n]


def slugify_worker(dask_worker):
    slug = re.sub("tcp://", "host_", dask_worker)
    slug = re.sub("\.", "_", slug)
    slug = re.sub(":", "-", slug)
    return slug


def match_geoms(ref_geom, geom_to_match, hydrogen=False):
    """
    See
        [1] 10.1021/ci400534h
        [2] 10.1021/acs.jcim.6b00516
    """

    ref_coords, _ = ref_geom.coords_by_type
    coords_to_match, inds_to_match = geom_to_match.coords_by_type
    atoms = ref_coords.keys()
    for atom in atoms:
        # Only match hydrogens if explicitly requested
        if atom == "H" and not hydrogen:
            continue
        print("atom", atom)
        ref_coords_for_atom = ref_coords[atom]
        coords_to_match_for_atom = coords_to_match[atom]
        # Pairwise distances between two collections
        # Atoms of ref_geom are along the rows, atoms of geom_to_match
        # along the columns.
        cd = cdist(ref_coords_for_atom, coords_to_match_for_atom)
        print(cd)
        # Hungarian method, row_inds are returned already sorted.
        row_inds, col_inds = linear_sum_assignment(cd)
        print("col_inds", col_inds)
        old_inds = inds_to_match[atom]
        new_inds = old_inds[col_inds]
        print("old_inds", old_inds)
        print("new_inds", new_inds)
        new_coords_for_atom = coords_to_match_for_atom[new_inds]
        # print(ref_coords_for_atom)
        # print(new_coords_for_atom)
        # print(ref_coords_for_atom-new_coords_for_atom)
        # Update coordinates
        print("old coords")
        c3d = geom_to_match.coords3d
        print(c3d)
        # Modify the coordinates directly
        c3d[old_inds] = new_coords_for_atom
        # coords_to_match[atom] = coords_to_match_for_atom[new_inds]


def check_for_stop_sign():
    stop_signs = ("stop", "STOP")
    stop_sign_found = False

    for ss in stop_signs:
        if os.path.exists(ss):
            print("Found stop sign. Stopping run.")
            os.remove(ss)
            stop_sign_found = True
    return stop_sign_found


def index_array_from_overlaps(overlaps, axis=1):
    """It is assumed that the overlaps between two points with indices
    i and j with (j > i) are computed and that i changes along the first
    axis (axis=0) and j changes along the second axis (axis=1).

    So the first row of the overlap matrix (overlaps[0]) should contain
    the overlaps between state 0 at index i and all states at index j.

    argmax along axis 1 returns the indices of the most overlapping states
    at index j with the states at index i, given by the item index in the
    indices array. E.g.:
        [0 1 3 2] indicates a root flip in a system with four states when
        going from index i to index j. Root 2 at i became root 3 at j and
        vice versa.
    """
    # indices = np.argmax(overlaps**2, axis=1)
    indices = np.argmax(np.abs(overlaps), axis=1)
    return indices


if __name__ == "__main__":
    print(load_geometry("hcn.xyz"))
    print(geoms_from_trj("cycle_040.trj"))
