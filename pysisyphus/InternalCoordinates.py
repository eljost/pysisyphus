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
    return bond_indices


def get_bending_indices(bond_indices):
    #for bi1 in bond_indices:
    bending_indices = list()
    for i in range(len(bond_indices)):
        #for bi2 in bond_indices:
        for j in range(i, len(bond_indices)):
            bi1 = bond_indices[i]
            bi2 = bond_indices[j]
            bond_set1 = set(bi1)
            bond_set2 = set(bi2)
            union = bond_set1 | bond_set2
            if len(union) == 3:
                central_set = bond_set1 & bond_set2
                assert len(central_set) == 1
                terminal1, terminal2 = union - central_set
                (central, ) = central_set
                bending_indices.append((terminal1, central, terminal2))
    logging.warning("No check for (nearly) linear angles!")
    return np.array(bending_indices, dtype=int)

if __name__ == "__main__":
    h2o_geom = geom_from_library("h2o.xyz") 
    h2o_inds = get_bond_indices(h2o_geom)
    assert len(h2o_inds) == 2
    h2o_bends = get_bending_indices(h2o_inds)
    print(h2o_bends)
    assert len(h2o_bends) == 1

    """
    benzene_geom = geom_from_library("benzene_bp86sto3g_opt.xyz")
    benzene_inds = get_bond_indices(benzene_geom)
    assert len(benzene_inds) == 12
    benzene_bends = get_bending_indices(benzene_inds)
    print(benzene_bends)
    """

    # See [2] for geometry
    fluorethylene = geom_from_library("fluorethylene.xyz")
    fe_inds = get_bond_indices(fluorethylene)
    assert len(fe_inds) == 5
    fe_bends = get_bending_indices(fe_inds)
    print(fe_bends+1, len(fe_bends))
    assert len(fe_bends) == 6

