#!/usr/bin/env python3

# [1] https://doi.org/10.1063/1.1515483
# [2] https://doi.org/10.1063/1.471864 delocalized internal coordinates

import itertools
import logging

import numpy as np
from scipy.spatial.distance import pdist

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
            term_ind = 1 - intersect_ind
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


def get_bond_B(geom, bond_ind):
    coords = geom.coords.reshape(-1, 3)
    n, m = bond_ind
    bond = coords[m] - coords[n]
    bond_length = np.linalg.norm(bond)
    bond_normed = bond / bond_length
    row = np.zeros_like(coords)
    # 1 / -1 correspond to the sign factor [1] Eq. 18
    row[m,:] = 1 * bond_normed
    row[n,:] = -1 * bond_normed
    row = row.flatten()
    return row


def get_angle_B(geom, angle_ind):
    def are_parallel(vec1, vec2, thresh=1e-6):
        rad = np.arccos(vec1.dot(vec2))
        return abs(rad) > (np.pi - thresh)
    coords = geom.coords.reshape(-1, 3)
    m, o, n = angle_ind
    u_dash = coords[m] - coords[o]
    v_dash = coords[n] - coords[o]
    u_norm = np.linalg.norm(u_dash)
    v_norm = np.linalg.norm(v_dash)
    u = u_dash / u_norm
    v = v_dash / v_norm
    rad = np.arccos(u.dot(v))
    # Eq. (24) in [1]
    if are_parallel(u, v):
        tmp_vec = np.array((1, -1, 1))
        par = are_parallel(u, tmp_vec) and are_parallel(v, tmp_vec)
        tmp_vec = np.array((-1, 1, 1)) if par else tmp_vec
        w_dash = np.cross(u, tmp_vec)
    else:
        w_dash = np.cross(u, v)
    w_norm = np.linalg.norm(w_dash)
    w = w_dash / w_norm
    uxw_dash = np.cross(u, w)
    wxv_dash = np.cross(w, v)
    uxw = uxw_dash / u_norm
    wxv = wxv_dash / v_norm

    row = np.zeros_like(coords)
    #                  |  m  |  n  |  o  |
    # -----------------------------------
    # sign_factor(amo) |  1  |  0  | -1 |
    # sign_factor(ano) |  0  |  1  | -1 |
    row[m,:] = uxw
    row[o,:] = -uxw - wxv
    row[n,:] = wxv
    row = row.flatten()
    return row


def get_dihedral_b(geom, dihedral_ind):
    print(dihedral_ind)
    coords = geom.coords.reshape(-1, 3)
    m, o, p, n = dihedral_ind
    u_dash = coords[m] - coords[o]
    v_dash = coords[n] - coords[p]
    w_dash = coords[p] - coords[o]
    u_norm = np.linalg.norm(u_dash)
    v_norm = np.linalg.norm(v_dash)
    w_norm = np.linalg.norm(w_dash)
    u = u_dash / u_norm
    v = v_dash / v_norm
    w = w_dash / w_norm
    phi_u = np.arccos(u.dot(w))
    phi_v = np.arccos(w.dot(v))
    uxw = np.cross(u, w)
    vxw = np.cross(v, w)
    dihed = np.arccos(uxw.dot(vxw)/(np.sin(phi_u)*np.sin(phi_v)))
    print(dihed)
    print(np.rad2deg(dihed))
    print()




if __name__ == "__main__":
    np.set_printoptions(suppress=True, precision=4)
    h2o_geom = geom_from_library("h2o.xyz") 
    h2o_inds = get_bond_indices(h2o_geom)
    assert len(h2o_inds) == 2
    h2o_bends = get_bending_indices(h2o_inds)
    assert len(h2o_bends) == 1
    [get_bond_B(h2o_geom, h2o_ind) for h2o_ind in h2o_inds]
    [get_angle_B(h2o_geom, h2o_bend) for h2o_bend in h2o_bends]

    benzene_geom = geom_from_library("benzene_bp86sto3g_opt.xyz")
    benzene_inds = get_bond_indices(benzene_geom)
    assert len(benzene_inds) == 12
    benzene_bends = get_bending_indices(benzene_inds)

    # See [2] for geometry
    fe_geom = geom_from_library("fluorethylene.xyz")
    fe_inds = get_bond_indices(fe_geom)
    assert len(fe_inds) == 5
    fe_bends = get_bending_indices(fe_inds)
    assert len(fe_bends) == 6
    fe_dihedrals = get_dihedral_indices(fe_inds, fe_bends)
    assert len(fe_dihedrals) == 4
    [get_dihedral_b(fe_geom, ind) for ind in fe_dihedrals]

    pt_geom = geom_from_library("h2o_pt.xyz")
    pt_inds = get_bond_indices(pt_geom)
    pt_bends = get_bending_indices(pt_inds)
    #print(pt_inds)
    #print(pt_bends)
    bbonds = [get_bond_B(pt_geom, ind) for ind in pt_inds]
    bbends = [get_angle_B(pt_geom, ind) for ind in pt_bends]
    bmat = np.array((*bbonds, *bbends))
    #print(bmat)
    #np.savetxt("pt", bmat)
