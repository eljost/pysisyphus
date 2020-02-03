import itertools as it
import numpy as np
from scipy.spatial.distance import pdist, squareform

from pysisyphus.elem_data import COVALENT_RADII as CR


def get_pair_covalent_radii(atoms):
    atoms = [a.lower() for a in atoms]
    cov_radii = np.array([CR[a] for a in atoms])
    pair_cov_radii = np.array([r1+r2 for r1, r2 in it.combinations(cov_radii, 2)])
    return pair_cov_radii


def get_bond_sets(atoms, coords3d, bond_factor=1.3):
    cdm = pdist(coords3d)
    # Generate indices corresponding to the atom pairs in the
    # condensed distance matrix cdm.
    atom_indices = list(it.combinations(range(len(coords3d)),2))
    atom_indices = np.array(atom_indices, dtype=int)
    cov_rad_sums = get_pair_covalent_radii(atoms)
    cov_rad_sums *= bond_factor
    bond_flags = cdm <= cov_rad_sums
    bond_indices = atom_indices[bond_flags]
    return bond_indices


def get_bond_mat(geom, bond_factor=1.3):
    cdm = pdist(geom.coords3d)
    pair_cov_radii = get_pair_covalent_radii(geom.atoms)
    bond_mat = squareform(cdm <= (pair_cov_radii * bond_factor))
    return bond_mat
