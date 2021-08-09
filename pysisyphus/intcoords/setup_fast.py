import itertools as it
import logging

import numpy as np
from sklearn.neighbors import KDTree

from pysisyphus.elem_data import COVALENT_RADII as CR
from pysisyphus.helpers_pure import log, timed
from pysisyphus.intcoords.valid import bend_valid, dihedral_valid


logger = logging.getLogger("internal_coords")

BOND_FACTOR = 1.3


def get_max_bond_dists(atoms, bond_factor, covalent_radii=None):
    if covalent_radii is None:
        cr = CR
    else:
        cr = {atom: covrad for atom, covrad in zip(atoms, covalent_radii)}

    unique_atoms = set(atoms)
    atom_pairs = it.combinations_with_replacement(unique_atoms, 2)
    max_bond_dists = {
        frozenset((from_, to_)): bond_factor * (cr[from_] + cr[to_])
        for from_, to_ in atom_pairs
    }
    return max_bond_dists


def find_bonds(atoms, coords3d, covalent_radii=None, bond_factor=BOND_FACTOR):
    atoms = [atom.lower() for atom in atoms]
    c3d = coords3d.reshape(-1, 3)
    if covalent_radii is None:
        covalent_radii = [CR[atom] for atom in atoms]
    cr = np.array(covalent_radii)

    max_bond_dists = get_max_bond_dists(atoms, bond_factor, covalent_radii=cr)
    radii = bond_factor * (cr.copy() + max(cr))
    kdt = KDTree(c3d)
    res, dists = kdt.query_radius(c3d, radii, return_distance=True)
    bonds_ = list()
    for i, (atom, bonds, dists_) in enumerate(zip(atoms, res, dists)):
        fsets = [frozenset((atom, atoms[a])) for a in bonds]
        ref_dists = np.array([max_bond_dists[fset] for fset in fsets])
        keep = np.logical_and(dists_ <= ref_dists, 0.1 <= dists_)
        for to_ in bonds[keep]:
            bonds_.append(frozenset((i, to_)))
    bonds = np.array([tuple((from_, to_)) for from_, to_ in set(bonds_)])
    return bonds


def find_bonds_for_geom(geom, bond_factor=BOND_FACTOR):
    return find_bonds(geom.atoms, geom.coords3d, geom.covalent_radii)


def get_bond_vec_getter(
    atoms, covalent_radii, bonds_for_inds, no_bonds_with=None, bond_factor=BOND_FACTOR
):
    max_bond_dists = get_max_bond_dists(
        atoms, bond_factor, covalent_radii=covalent_radii
    )
    max_bond_dists_for_inds = [
        np.array([max_bond_dists[frozenset((atoms[ind], atom_))] for atom_ in atoms])
        for ind in bonds_for_inds
    ]
    if no_bonds_with is None:
        # List of empty lists
        no_bonds_with = [[] * len(bonds_for_inds)]

    # Set it to a negative value, so the calculated distance, which is always positive
    # can't be smaller than this value.
    for mbd, nbw in zip(max_bond_dists_for_inds, no_bonds_with):
        mbd[nbw] = -1

    all_inds = np.arange(len(atoms))

    def get_bond_vecs(coords, return_bonded_inds=False):
        coords3d = coords.reshape(-1, 3)
        all_bond_vecs = list()
        all_bonded_inds = list()
        for ind, max_dists in zip(bonds_for_inds, max_bond_dists_for_inds):
            distance_vecs = coords3d - coords3d[ind]
            distances = np.linalg.norm(distance_vecs, axis=1)
            # Set 0.0 distance of atom with itself to a high value to not form
            # and ind-ind bond.
            distances[ind] = 10_000
            bond_mask = distances <= max_dists
            all_bond_vecs.append(distance_vecs[bond_mask])
            all_bonded_inds.append(all_inds[bond_mask])

        if return_bonded_inds:
            return all_bond_vecs, all_bonded_inds
        else:
            return all_bond_vecs

    return get_bond_vecs


def find_bends(coords3d, bonds, min_deg, max_deg, logger=None):
    bond_dict = {}
    bonds = [tuple(bond) for bond in bonds]
    for from_, to_ in bonds:
        bond_dict.setdefault(from_, list()).append(to_)
        bond_dict.setdefault(to_, list()).append(from_)

    bend_set = set()
    for bond in bonds:
        from_, to_ = bond
        from_neighs = set(bond_dict[from_]) - set((to_,))
        to_neighs = set(bond_dict[to_]) - set((from_,))
        bend_candidates = [(neigh,) + bond for neigh in from_neighs] + [
            bond + (neigh,) for neigh in to_neighs
        ]
        for indices in bend_candidates:
            if (
                (not bend_valid(coords3d, indices, min_deg, max_deg))
                or (indices in bend_set)
                or (indices[::-1] in bend_set)
            ):
                continue
            bend_set.add(indices)
    return [list(bend) for bend in bend_set]


def find_dihedrals(coords3d, bonds, bends, max_deg, logger=None):
    bond_dict = {}
    bonds = [tuple(bond) for bond in bonds]
    for from_, to_ in bonds:
        bond_dict.setdefault(from_, list()).append(to_)
        bond_dict.setdefault(to_, list()).append(from_)

    dihedral_set = set()
    for bend in bends:
        bend = tuple(bend)
        from_, central, to_ = bend
        from_neighs = set(bond_dict[from_]) - set((to_, central))
        to_neighs = set(bond_dict[to_]) - set((from_, central))
        dihedral_candidates = [(neigh,) + bend for neigh in from_neighs] + [
            bend + (neigh,) for neigh in to_neighs
        ]
        for indices in dihedral_candidates:
            if not dihedral_valid(coords3d, indices, deg_thresh=max_deg):
                continue
            dihedral_set.add(indices)
    return [list(dihedral) for dihedral in dihedral_set]


def find_bonds_bends(geom, bond_factor=BOND_FACTOR, min_deg=15, max_deg=175):
    log(logger, "Starting detection of bonds and bends.")
    bonds = find_bonds_for_geom(geom, bond_factor=bond_factor)
    log(logger, f"Found {len(bonds)} bonds.")
    bends = find_bends(geom.coords3d, bonds, min_deg=min_deg, max_deg=max_deg)
    log(logger, f"Found {len(bends)} bends.")
    return bonds, bends


@timed(logger)
def find_bonds_bends_dihedrals(geom, bond_factor=BOND_FACTOR, min_deg=15, max_deg=175):
    log(logger, f"Detecting bonds, bends and dihedrals for {len(geom.atoms)} atoms.")
    bonds, bends = find_bonds_bends(
        geom, bond_factor=bond_factor, min_deg=min_deg, max_deg=max_deg
    )
    proper_dihedrals = find_dihedrals(geom.coords3d, bonds, bends, max_deg)
    log(logger, f"Found {len(proper_dihedrals)} proper dihedrals.")
    return bonds, bends, proper_dihedrals
