import itertools as it
import numpy as np
from scipy.spatial.distance import pdist

from pysisyphus.elem_data import COVALENT_RADII as CR


def get_cov_radii_sum_array(atoms, coords):
    coords3d = coords.reshape(-1, 3)
    atom_indices = list(it.combinations(range(len(coords3d)),2))
    atom_indices = np.array(atom_indices, dtype=int)
    cov_rad_sums = list()
    for i, j in atom_indices:
        atom1 = atoms[i].lower()
        atom2 = atoms[j].lower()
        cov_rad_sums.append(CR[atom1] + CR[atom2])
    cov_rad_sums = np.array(cov_rad_sums)
    return cov_rad_sums


def get_bond_sets(atoms, coords3d, bond_factor=1.3):
    cdm = pdist(coords3d)
    # Generate indices corresponding to the atom pairs in the
    # condensed distance matrix cdm.
    atom_indices = list(it.combinations(range(len(coords3d)),2))
    atom_indices = np.array(atom_indices, dtype=int)
    cov_rad_sums = get_cov_radii_sum_array(atoms, coords3d.flatten())
    cov_rad_sums *= bond_factor
    bond_flags = cdm <= cov_rad_sums
    bond_indices = atom_indices[bond_flags]
    return bond_indices
