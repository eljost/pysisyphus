#!/usr/bin/env python3

import argparse
import logging
import itertools as it
import re
import sys

import numpy as np

from pysisyphus.helpers import geoms_from_trj, geom_from_xyz_file
from pysisyphus.InternalCoordinates import RedundantCoords


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("trj")
    parser.add_argument("filter")
    parser.add_argument("--first", type=int, metavar="N", default=None,
                        help="Only consider first N geometries in .trj."
    )
    parser.add_argument("--cov_rad_factor", type=float, default=1.3,
                        help="Factor to determine bonds."
    )

    return parser.parse_args(args)


def parse_filter(raw_filter):
    split = raw_filter.strip().split()
    mobjs = [re.match("(?P<not_present>!?)\((?P<atoms>[a-zA-Z-]+)\)", s) for s in split]
    filters = list()
    for m in mobjs:
        internal = m["atoms"]
        is_present = not bool(m["not_present"])
        # internal = set(internal.upper().split("-"))
        internal = tuple(sorted(internal.upper().split("-")))
        filters.append((is_present, internal))
    return filters


def get_unique_internals(geom):
    attrs = ("bond_indices", "bending_indices", "dihedral_indices")
    atoms_arr = np.array(geom.atoms)
    unique_internals = list()
    for attr in attrs:
        indices = getattr(geom.internal, attr)
        # unique = set([tuple(atoms_arr[inds]) for inds in indices])
        unique = set([tuple(sorted(atoms_arr[inds])) for inds in indices])
        unique_internals.extend(unique)
    return unique_internals


def run():
    args = parse_args(sys.argv[1:])

    internal_kwargs = {
        "cov_rad_factor": args.cov_rad_factor,
    }
    geoms = geoms_from_trj(args.trj, first=args.first, coord_type="redund",
                           internal_kwargs=internal_kwargs)
    atoms = geoms[0].atoms
    filters = parse_filter(args.filter)
    print("Filters", filters)
    # print(ai)
    # geoms = [geoms[19], ]
    gui = get_unique_internals
    for i, geom in enumerate(geoms):
        internals = get_unique_internals(geom)
        valid = all([is_present == (internal in internals)
                     for is_present, internal in filters]
        )
        # print(valid)
        if not valid:
            print(f"geom {i+1} is not valid.")
            print(internals)
    # import pdb; pdb.set_trace()
