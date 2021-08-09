import argparse
import re
import sys

import numpy as np

from pysisyphus.helpers import geoms_from_trj
from pysisyphus.xyzloader import write_geoms_to_trj


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("trj")
    parser.add_argument("filter")
    parser.add_argument("--first", type=int, metavar="N", default=None,
                        help="Only consider first N geometries in .trj."
    )

    return parser.parse_args(args)


def parse_filter(raw_filter):
    split = raw_filter.strip().split()
    mobjs = [re.match("(?P<not_present>!?)\((?P<atoms>[a-zA-Z-]+)\)", s) for s in split]
    filters = list()
    for m in mobjs:
        internal = m["atoms"]
        is_present = not bool(m["not_present"])
        internal = tuple(sorted(internal.upper().split("-")))
        filters.append((is_present, internal))
    return filters


def get_unique_internals(geom):
    attrs = ("bond_atom_indices", "bend_atom_indices", "dihedral_atom_indices")
    atoms_arr = np.array(geom.atoms)
    unique_internals = list()
    for attr in attrs:
        indices = getattr(geom.internal, attr)
        unique = set([tuple(sorted(atoms_arr[inds])) for inds in indices])
        unique_internals.extend(unique)
    return unique_internals


def run():
    args = parse_args(sys.argv[1:])
    geoms = geoms_from_trj(args.trj, first=args.first, coord_type="redund",)
    filters = parse_filter(args.filter)
    print("Filters", filters)

    valid_geoms = list()
    invalid_geoms = list()
    for i, geom in enumerate(geoms):
        internals = get_unique_internals(geom)
        valid = all([is_present == (internal in internals)
                     for is_present, internal in filters]
        )
        if not valid:
            print(f"geom {i+1} is not valid.")
            invalid_geoms.append(geom)
        else:
            valid_geoms.append(geom)
    print()
    print(f"Found {len(valid_geoms)} valid geometries.")
    print(f"Found {len(invalid_geoms)} invalid geometries.")

    write_geoms_to_trj(valid_geoms, "filtered_valid.trj")
    write_geoms_to_trj(invalid_geoms, "filtered_invalid.trj")
