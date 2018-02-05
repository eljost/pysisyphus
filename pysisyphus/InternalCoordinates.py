#!/usr/bin/env python3

# [1] https://doi.org/10.1063/1.1515483
# [2] https://doi.org/10.1063/1.471864 delocalized internal coordinates

import itertools
import logging

import numpy as np
from scipy.spatial.distance import pdist, squareform

from pysisyphus.helpers import geom_from_library
from pysisyphus.elem_data import COVALENT_RADII as CR

def get_bond_indices(geom, factor=1.3):
    """
    Default factor taken from [1] A.1.
    """
    coords = geom.coords.reshape(-1, 3)
    # Condensed distance matrix
    cdm = pdist(coords)
    # Generate indices corresponding to the atom pairs in
    # condensed distance matrix cdm.
    atom_indices = list(itertools.combinations(range(len(coords)),2))
    atom_indices = np.array(atom_indices, dtype=int)
    cov_rad_sums = list()
    for i, j in atom_indices:
        atom1 = geom.atoms[i].lower()
        atom2 = geom.atoms[j].lower()
        cov_rad1 = CR[atom1]
        cov_rad2 = CR[atom2]
        cov_rad_sum = factor * (cov_rad1 + cov_rad2)
        cov_rad_sums.append(cov_rad_sum)
    cov_rad_sums = np.array(cov_rad_sums)
    bond_flags = cdm <= cov_rad_sums
    bond_indices = atom_indices[bond_flags]
    logging.warning("No check for hydrogen bonds or disconnected fragments!")
    # Set of sets with frozenset?
    return bond_indices


def sort_by_central(set1, set2):
    """Determines a common index in two sets and returns a length 3
    tuple with the central index at the middle position and the two
    terminal indices as first and last indices."""
    central_set = set1 & set2
    union = set1 | set2
    assert len(central_set) == 1
    terminal1, terminal2 = union - central_set
    (central, ) = central_set
    return (terminal1, central, terminal2), central


def get_bending_indices(bond_indices):
    bond_sets = {frozenset(bi) for bi in bond_indices}
    bending_indices = list()
    for bond_set1, bond_set2 in itertools.combinations(bond_sets, 2):
        union = bond_set1 | bond_set2
        if len(union) == 3:
            as_tpl, _ = sort_by_central(bond_set1, bond_set2)
            bending_indices.append(as_tpl)
    logging.warning("No check for (nearly) linear angles!"
                    "No additional orthogonal bending coordinates.")
    return np.array(bending_indices, dtype=int)


def get_dihedral_indices(bond_inds, bend_inds):
    dihedral_inds = list()
    dihedral_sets = list()
    for bond, bend in itertools.product(bond_inds, bend_inds):
        central = bend[1]
        bes = set((bend[0], bend[2]))
        bois = set(bond)
        if (len(bes & bois) == 1) and (central not in bois):
            (intersect,)  = set(bond) & set(bend)
            intersect_ind = list(bond).index(intersect)
            term_ind = 0 if intersect_ind else 1
            terminal = bond[term_ind]
            if intersect == bend[0]:
                dihedral_ind = [terminal] + list(bend)
            else:
                dihedral_ind = list(bend) + [terminal]
            dihedral_set = set(dihedral_ind)
            if dihedral_set not in dihedral_sets:
                dihedral_inds.append(dihedral_ind)
                dihedral_sets.append(dihedral_set)
    logging.warning("No check for dihedrals near 180° or -180°!")
    return np.array(dihedral_inds)



if __name__ == "__main__":
    h2o_geom = geom_from_library("h2o.xyz") 
    h2o_inds = get_bond_indices(h2o_geom)
    assert len(h2o_inds) == 2
    h2o_bends = get_bending_indices(h2o_inds)
    assert len(h2o_bends) == 1

    benzene_geom = geom_from_library("benzene_bp86sto3g_opt.xyz")
    benzene_inds = get_bond_indices(benzene_geom)
    assert len(benzene_inds) == 12
    benzene_bends = get_bending_indices(benzene_inds)

    # See [2] for geometry
    fluorethylene = geom_from_library("fluorethylene.xyz")
    fe_inds = get_bond_indices(fluorethylene)
    assert len(fe_inds) == 5
    fe_bends = get_bending_indices(fe_inds)
    assert len(fe_bends) == 6
    fe_dihedrals = get_dihedral_indices(fe_inds, fe_bends)
    print(fe_dihedrals)
