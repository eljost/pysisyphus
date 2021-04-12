#!/usr/bin/env python

import argparse
import sys

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_loader
from pysisyphus.helpers_pure import full_expand
from pysisyphus.intcoords.setup import get_fragments, get_bond_sets


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("fn")
    parser.add_argument("bonds", type=int, nargs="+")

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
    args = parse_args(sys.argv[1:])

    geom = geom_loader(args.fn)
    rm_bonds_ = args.bonds
    assert (len(rm_bonds_) > 0) and (len(rm_bonds_) % 2) == 0
    rm_bonds = list()
    for i in range(len(rm_bonds_) // 2):
        rm_bond = rm_bonds_[2 * i : 2 * i + 2]
        rm_bond.sort()
        rm_bonds.append(rm_bond)
    print(rm_bonds)

    bonds = get_bond_sets(geom.atoms, geom.coords3d).tolist()
    for rm_bond in rm_bonds:
        try:
            bonds.remove(rm_bond)
        except ValueError:
            print(f"Bond {rm_bond} not in detected bonds!")

    fragments = get_fragments(geom.atoms, geom.coords, bond_inds=bonds)
    print(f"Found {len(fragments)} fragments")

    atoms = geom.atoms
    dummy_atoms = ["X"] * len(atoms)
    coords3d = geom.coords3d
    xyzs = list()
    for i, frag in enumerate(fragments):
        comped = compress(frag)
        print(f"{len(frag)} atoms, compressed:\n\t{comped}")
        # frag_atoms = [atoms[i] for i in frag]
        frag_atoms = dummy_atoms.copy()
        for j in frag:
            frag_atoms[j] = atoms[j]
        frag_geom = Geometry(frag_atoms, coords3d)
        fn = f"frag_{i:03d}.xyz"
        as_xyz = frag_geom.as_xyz()
        with open(fn, "w") as handle:
            handle.write(as_xyz)
        xyzs.append(as_xyz)
    with open("fragments.trj", "w") as handle:
        handle.write("\n".join(xyzs))


if __name__ == "__main__":
    run()
    # res = compress((10, 11, 12, 13, 17, 18, 19, 26))
    # print(res)
