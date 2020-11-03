import itertools as it

import numpy as np
from sklearn.neighbors import KDTree

from pysisyphus.elem_data import COVALENT_RADII as CR
from pysisyphus.intcoords.valid import bend_valid


def find_bonds(geom, bond_factor=1.3):
    atoms = [atom.lower() for atom in geom.atoms]
    c3d = geom.coords3d
    unique_atoms = set(atoms)
    atom_pairs = it.combinations_with_replacement(unique_atoms, 2)
    max_bond_dists = {
        frozenset((from_, to_)): bond_factor * (CR[from_] + CR[to_])
        for from_, to_ in atom_pairs
    }

    cr = geom.covalent_radii
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


def find_bonds_bends(geom, bond_factor=1.3):
    bonds = find_bonds(geom, bond_factor=bond_factor)
    bends = find_bends(geom.coords3d, bonds, min_deg=15, max_deg=175)
    return bonds, bends
