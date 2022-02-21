#!/usr/bin/env python

import argparse
import itertools as it
import sys

from pysisyphus.elem_data import ATOMIC_NUMBERS
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_loader
from pysisyphus.helpers_pure import full_expand, highlight_text
from pysisyphus.intcoords.setup import get_fragments, get_bond_sets


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("fn")
    break_group = parser.add_mutually_exclusive_group(required=True)
    break_group.add_argument(
        "--bonds",
        type=int,
        nargs="+",
        help="List of pairs of atom indices. Each pair is a bond that "
        "is broken, to generate the fragments.",
    )
    break_group.add_argument(
        "--bonds-from",
        type=int,
        nargs="+",
        help="List of atom indices. Break all bonds involving any of this atoms.",
    )
    break_group.add_argument(
        "--tmc",
        action="store_true",
        help="Break all bonds involving transition metal atoms.",
    )
    return parser.parse_args(args)


def compress(inds, check=True):
    inds = list(inds)
    inds.sort()

    start, *rest = inds
    compressed = str(start)

    prev_ind = start
    cur_range = [
        start,
    ]
    for i, ind in enumerate(rest):
        if ind != prev_ind + 1:
            if len(cur_range) > 1:
                compressed += f"..{cur_range[-1]+1},{ind}"
            else:
                compressed += f",{ind}"
            cur_range = [
                ind,
            ]
        else:
            cur_range.append(ind)
        prev_ind = ind
    else:
        if len(cur_range) > 1:
            compressed += f"..{cur_range[-1]+1}"
    if check:
        expd = full_expand(compressed)
        assert len(expd) == len(inds)
    return compressed


def run():
    print(highlight_text("frag.py") + "\n")
    args = parse_args(sys.argv[1:])

    geom = geom_loader(args.fn)

    # Determine bond topology
    bonds = get_bond_sets(geom.atoms, geom.coords3d).tolist()

    rm_bonds_from = set()
    if rm_bonds_ := args.bonds:
        # At least one bond must be present and the atom indices must be given in pairs
        assert (len(rm_bonds_) > 0) and (len(rm_bonds_) % 2) == 0
        rm_bonds = list()
        for i in range(len(rm_bonds_) // 2):
            rm_bond = rm_bonds_[2 * i : 2 * i + 2]
            rm_bond.sort()
            rm_bonds.append(rm_bond)
    elif args.bonds_from:
        rm_bonds_from = set(args.bonds_from)
    elif args.tmc:
        print("Delete all bonds involving transition metals")
        trans_metal_nums = list(it.chain(
            #       3d             4d             5d             6d
            *(range(21, 31), range(39, 49), range(57, 81), range(89, 113))
        ))

        def is_trans_metal(atomic_num):
            return atomic_num in trans_metal_nums

        atom_nums = [ATOMIC_NUMBERS[atom.lower()] for atom in geom.atoms]
        trans_metal_atoms = [
            i for i, atom_num in enumerate(atom_nums) if is_trans_metal(atom_num)
        ]
        rm_bonds_from = set(trans_metal_atoms)
    else:
        raise Exception("How did I get here?")

    if rm_bonds_from:
        rm_atoms = [geom.atoms[i] for i in rm_bonds_from]
        rm_str = ", ".join([f"{i}{atom}" for i, atom in zip(rm_bonds_from, rm_atoms)])
        print(f"Delete all bonds involving atoms: {rm_str}")
        rm_bonds = [bond for bond in bonds if rm_bonds_from & set(bond)]
    print()

    assert rm_bonds

    print(f"{len(rm_bonds)} bond(s) to be deleted:")
    for i, bond in enumerate(rm_bonds):
        print(f"\t{i:02d}: {bond}")

    for rm_bond in rm_bonds:
        try:
            bonds.remove(rm_bond)
        except ValueError:
            print(f"Bond {rm_bond} not in detected bonds!")
    print()

    fragments = get_fragments(geom.atoms, geom.coords, bond_inds=bonds)
    single_atom_frags = set(range(len(geom.atoms))) - set(it.chain(*fragments))
    if single_atom_frags:
        print(f"Found {len(single_atom_frags)} single atom fragments.")
    fragments += [
        [
            single_atom,
        ]
        for single_atom in single_atom_frags
    ]
    print(f"Found {len(fragments)} fragments")

    atoms = geom.atoms
    dummy_atoms = ["X"] * len(atoms)
    coords3d = geom.coords3d
    xyzs = list()
    frags_compressed = list()
    for i, frag in enumerate(fragments):
        compressed = compress(frag, check=True)
        frags_compressed.append(compressed)
        sop = "s" if len(frag) > 1 else ""
        print(f"{len(frag)} atom{sop}, compressed:\n\t{compressed}")
        frag_atoms = dummy_atoms.copy()
        for j in frag:
            frag_atoms[j] = atoms[j]
        frag_geom = Geometry(frag_atoms, coords3d)
        fn = f"frag_{i:03d}.xyz"
        as_xyz = frag_geom.as_xyz()
        with open(fn, "w") as handle:
            handle.write(as_xyz)
        xyzs.append(as_xyz)

    # Do a round trip to see if all fragments add up to the total number of atoms
    expanded = set(full_expand(",".join(frags_compressed)))
    assert len(expanded) == len(geom.atoms)

    with open("fragments.trj", "w") as handle:
        handle.write("\n".join(xyzs))


if __name__ == "__main__":
    run()
