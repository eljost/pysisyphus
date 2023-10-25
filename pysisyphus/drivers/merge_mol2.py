import argparse
import copy
import sys
from typing import Optional
import warnings

import numpy as np

from pysisyphus.constants import BOHR2ANG
from pysisyphus.drivers.merge import merge_with_frozen_geom
from pysisyphus.helpers import geom_loader
from pysisyphus.io.mol2 import parse_mol2, dict_to_mol2_string


def delete_atoms_bonds_inplace(
    as_dict: dict, inds: list[int], atom_offset: int = 0, bond_offset: int = 0
) -> dict[int, int]:
    """Update atom_ids and bonds in mol2-dict inplace.

    Parameter
    ---------
    as_dict
        mol2-dict as returned from pysisyphus.io.mol2.parse_mol2.
    inds
        List of positive integer atom_ids to be deleted.
    atom_offset
        Integer >= 0. atom_ids will be shifted by this number.
    bond_offset
        Integer >= 0. bond_ids will be shifted by this number.

    Returns
    -------
    atom_map
        Dictionary w/ origin atom_ids as keys and updated atom_ids as values.
    """
    to_del = set(inds)

    # Atoms
    axs = as_dict["atoms_xyzs"]
    axs_mod = list()
    atoms_del = 0
    atom_map = {}
    for ax in axs:
        if ax["atom_id"] in to_del:
            warnings.warn("Charges were not updated after deleting atoms!")
            print("Deleted", ax)
            atoms_del += 1
            continue
        atom_id = ax["atom_id"]
        mod_atom_id = atom_id - atoms_del + atom_offset
        # Store updated atom_id in dict with original atom_id as key,
        # so we can later acces them to update the bond atom_ids.
        atom_map[atom_id] = mod_atom_id
        ax["atom_id"] = mod_atom_id
        axs_mod.append(ax)
    as_dict["atoms_xyzs"] = axs_mod

    # Bonds
    bonds_mod = list()
    bonds_del = 0
    for bond in as_dict["bond"]:
        bond_inds = set((bond["origin_atom_id"], bond["target_atom_id"]))
        if bond_inds & to_del:
            print("Deleted", bond)
            bonds_del += 1
            continue
        bond["bond_id"] -= bonds_del
        bond["bond_id"] += bond_offset
        bond["origin_atom_id"] = atom_map[bond["origin_atom_id"]]
        bond["target_atom_id"] = atom_map[bond["target_atom_id"]]
        bonds_mod.append(bond)
    as_dict["bond"] = bonds_mod
    return atom_map


def update_coords_inplace(as_dict: dict, coords3d: np.ndarray):
    """Update xyz coordinates in mol2 dict inplace."""
    axs = as_dict["atoms_xyzs"]
    assert len(axs) == len(coords3d)

    for ax, nc3d in zip(axs, coords3d):
        ax["xyz"] = (nc3d * BOHR2ANG).tolist()


def merge_mol2_dicts(as_dict1: dict, as_dict2: dict) -> dict:
    """Merge two mol2 dict. Atoms, bonds and associated numbers will be updated.

    Parameters
    ----------
    as_dict1
        Mol2 dict.
    as_dict2
        Mo2 dict.

    Returns
    -------
    merged
        Merged mol2 dict. Most entries will be copied from as_dict1 but atoms,
        bonds and associated counts will be updated with entries from as_dict2.
    """
    # Atoms
    atoms_xyzs = as_dict1["atoms_xyzs"] + as_dict2["atoms_xyzs"]

    # Bonds
    bond = as_dict1["bond"] + as_dict2["bond"]

    merged = copy.deepcopy(as_dict1)
    merged.update(
        {
            "atoms_xyzs": atoms_xyzs,
            "num_atoms": len(atoms_xyzs),
            "bond": bond,
            "num_bonds": len(bond),
            "mol_name": "merged",
        }
    )
    return merged


def merge_mol2(
    fn1: str,
    fn2: str,
    del1: Optional[list[int]] = None,
    del2: Optional[list[int]] = None,
    bonds: Optional[list[int]] = None,
    new_coords: Optional[np.ndarray] = None,
) -> dict:
    if del1 is None:
        del1 = list()
    if del2 is None:
        del2 = list()
    if bonds is None:
        bonds = list()
    bonds = np.array(bonds, dtype=int).reshape(-1, 2)

    dict1 = parse_mol2(fn1).as_dict()
    # dict1 = m1.as_dict()
    dict2 = parse_mol2(fn2).as_dict()
    # dict2 = m2.as_dict()

    # Delete atoms and assoicated bonds
    atom_map1 = delete_atoms_bonds_inplace(dict1, del1)
    print()
    natoms1 = len(dict1["atoms_xyzs"])
    nbonds1 = len(dict1["bond"])

    atom_map2 = delete_atoms_bonds_inplace(
        dict2, del2, atom_offset=natoms1, bond_offset=nbonds1
    )

    if new_coords is not None:
        coords1 = new_coords[:natoms1]
        coords2 = new_coords[natoms1:]
        assert len(coords1) + len(coords2) == len(new_coords)
        update_coords_inplace(dict1, coords1)
        update_coords_inplace(dict2, coords2)
        print("Updated coordinates")

    # Add new bond(s) to dict2
    bond2 = dict2["bond"]
    nbonds2 = len(bond2)
    for from_, to_ in bonds:
        # Use updated atom_ids for the bond target/origin
        origin = atom_map1[from_]
        target = atom_map2[to_]
        nbonds2 += 1
        bond_type = 1
        new_bond = {
            "bond_id": nbonds1 + nbonds2,
            "bond_type": bond_type,
            "origin_atom_id": origin,
            "target_atom_id": target,
        }
        print(f"Added new bond '{new_bond}'")
        warnings.warn(f"Set bond_type to '{bond_type}'.")
        bond2.append(new_bond)

    merged = merge_mol2_dicts(dict1, dict2)
    return merged


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="Merge two mol2 files w/ atom deletion & bond formation."
    )

    parser.add_argument("fn1", help="Name of frozen mol2 file.")
    parser.add_argument("fn2", help="Name of mobile mol2 file.")
    parser.add_argument(
        "--del1",
        nargs="+",
        type=int,
        help="1-based atom ids to be deleted from fn1.",
    )
    parser.add_argument(
        "--del2",
        nargs="+",
        type=int,
        help="1-based atom ids to be deleted from fn2.",
    )
    parser.add_argument(
        "--bonds",
        nargs="+",
        type=int,
        help="1-based atom id pairs (id1, id2), between which bonds are formed."
        "The original atom ids must be used, regardless of any atom deletions.",
    )
    parser.add_argument(
        "--out", default="merged.mol2", help="Name of the final mol2 file."
    )

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    fn1 = args.fn1
    fn2 = args.fn2
    del1 = list(sorted(args.del1))
    del2 = list(sorted(args.del2))
    bonds = args.bonds
    out = args.out

    # Make 0-based indices from 1-based mol2 indices
    del10 = [d - 1 for d in args.del1]
    del20 = [d - 1 for d in args.del2]
    bonds0 = np.array(bonds, dtype=int).reshape(-1, 2) - 1

    geom1 = geom_loader(fn1)
    geom2 = geom_loader(fn2)
    # Get new coordinates by merging both fragments.
    # Function takes 0-based indices
    new_geom, *_ = merge_with_frozen_geom(geom1, geom2, bonds0, del10, del20)
    new_coords3d = new_geom.coords3d

    # Merge mol2 files w/ updated coordinates
    # Takes 1-based indices, bad?!
    merged = merge_mol2(fn1, fn2, del1, del2, bonds, new_coords3d)

    # Render mol2 string from merged data
    rendered = dict_to_mol2_string(merged)

    # Dump merged mol2 dict to file
    with open(out, "w") as handle:
        handle.write(rendered)
    print(f"Dumped merged geometry to '{out}'.")


if __name__ == "__main__":
    run()
