import os
from pathlib import Path

import numpy as np
import scipy as sp

from pysisyphus.constants import ANG2BOHR
from pysisyphus.Geometry import Geometry
from pysisyphus.xyzloader import parse_xyz_file, parse_trj_file


# Taken from periodictable-1.5.0
MASS_DICT = {
    'n': 14.0067, 'h': 1.00794, 'he': 4.002602, 'li': 6.941, 'be': 9.012182,
    'b': 10.811, 'c': 12.0107, 'o': 15.9994, 'f': 18.9984032, 'ne': 20.1797,
    'na': 22.98977, 'mg': 24.305, 'al': 26.981538, 'si': 28.0855, 'p': 30.973761,
    's': 32.065, 'cl': 35.453, 'ar': 39.948, 'k': 39.0983, 'ca': 40.078,
    'sc': 44.95591, 'ti': 47.867, 'v': 50.9415, 'cr': 51.9961, 'mn': 54.938049,
    'fe': 55.845, 'co': 58.9332, 'ni': 58.6934, 'cu':63.546, 'zn': 65.409,
    'ga': 69.723, 'ge': 72.64, 'as': 74.9216, 'se': 78.96, 'br': 79.904,
    'kr': 83.798,'rb': 85.4678, 'sr': 87.62, 'y': 88.90585, 'zr': 91.224,
    'nb': 92.90638, 'mo': 95.94, 'tc': 98, 'ru': 101.07, 'rh': 102.9055,
    'pd': 106.42, 'ag': 107.8682, 'cd': 112.411, 'in': 114.818, 'sn': 118.71,
    'sb': 121.76, 'te': 127.6, 'i': 126.90447, 'xe': 131.293, 'cs': 132.90545,
    'ba': 137.327, 'la': 138.9055, 'ce': 140.116, 'pr': 140.90765, 'nd': 144.24,
    'pm': 145, 'sm': 150.36, 'eu': 151.964, 'gd': 157.25, 'tb': 158.92534,
    'dy': 162.5, 'ho': 164.93032, 'er': 167.259, 'tm': 168.93421, 'yb': 173.04,
    'lu': 174.967, 'hf': 178.49, 'ta': 180.9479, 'w': 183.84, 're': 186.207,
    'os': 190.23, 'ir': 192.217, 'pt': 195.078, 'au': 196.96655, 'hg': 200.59,
    'tl': 204.3833, 'pb': 207.2, 'bi': 208.98038, 'po': 209, 'at': 210, 'rn': 222,
    'fr': 223, 'ra': 226, 'ac': 227, 'th': 232.0381, 'pa': 231.03588, 'u': 238.02891,
    'np': 237, 'pu': 244, 'am': 243, 'cm': 247, 'bk': 247, 'cf': 251, 'es': 252,
    'fm': 257, 'md': 258, 'no': 259, 'lr': 262, 'rf': 261, 'db': 262, 'sg': 266,
    'bh': 264, 'hs': 277, 'mt': 268, 'ds': 281, 'rg': 272, 'cn': 285, 'nh': 286,
    'fl': 289,'mc': 289, 'lv': 293, 'ts': 294, 'og': 294
}


def geom_from_xyz_file(xyz_fn):
    atoms, coords = parse_xyz_file(xyz_fn)
    coords *= ANG2BOHR
    geom = Geometry(atoms, coords.flatten())
    return geom


def geom_from_library(xyz_fn):
    this_dir = os.path.abspath(os.path.dirname(__file__))
    xyz_dir = Path(this_dir) / "../xyz_files/"
    atoms, coords = parse_xyz_file(xyz_dir / xyz_fn)
    coords *= ANG2BOHR
    geom = Geometry(atoms, coords.flatten())
    return geom


def geoms_from_trj(trj_fn):
    geoms = [Geometry(atoms, coords.flatten()*ANG2BOHR)
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
