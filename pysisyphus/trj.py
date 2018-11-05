#!/usr/bin/env python3

import argparse
import copy
import itertools
import os
from pathlib import Path
from pprint import pprint
import re
import sys

from natsort import natsorted
import yaml

from pysisyphus.cos import *
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_xyz_file, geoms_from_trj, procrustes
from pysisyphus.calculators.IDPP import idpp_interpolate
from pysisyphus.xyzloader import write_geoms_to_trj


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("fns", nargs="+",
            help="Filenames of .xyz and/or .trj files (xyz and trj can be mixed)."
    )

    action_group = parser.add_mutually_exclusive_group(required=True)
    action_group.add_argument("--between", type=int,
                    help="Interpolate additional images."
    )
    action_group.add_argument("--align", action="store_true",
                    help="Align geometries onto the first geometry read from "
                         " multiple .xyz  or one .trj file."
    )
    action_group.add_argument("--split", action="store_true",
                    help="Split a supplied .trj file in multiple .xyz files."
    )
    action_group.add_argument("--reverse", action="store_true",
                    help="Reverse a .trj file."
    )
    action_group.add_argument("--cleantrj", action="store_true",
                    help="Clean a .trj file."
    )
    action_group.add_argument("--spline", action="store_true",
                    help="Evenly redistribute geometries along a splined path."
    )
    action_group.add_argument("--first", type=int,
                    help="Copy the first N geometries to a new file."
    )
    action_group.add_argument("--every", type=int,
                    help="Create new .trj with every N-th geometry. "
                         "Always include the first and last point."
    )
    action_group.add_argument("--bohr2ang", action="store_true",
                    help="Convert geometries from Bohr to Angstrom."
    )

    parser.add_argument("--idpp", action="store_true",
        help="Use Image Dependent Pair Potential instead "
             "of simple linear interpolation.")
    return parser.parse_args()


def read_geoms(xyz_fns, comments=False):
    if isinstance(xyz_fns, str):
        xyz_fns = [xyz_fns, ]

    geoms = list()
    for fn in xyz_fns:
        if fn.endswith(".xyz"):
            geom = [geom_from_xyz_file(fn), ]
        elif fn.endswith(".trj"):
            geom = geoms_from_trj(fn)
        else:
            raise Exception("Only .xyz and .trj files are supported!")
        geoms.extend(geom)
    return geoms


def get_geoms(xyz_fns, idpp=False, between=0, comments=False):
    """Returns a list of Geometry objects."""

    geoms = read_geoms(xyz_fns, comments=comments)

    print(f"Read {len(geoms)} geometries.")

    # Do IDPP interpolation if requested,
    trj = ""
    xyz_per_image = list()
    if idpp:
        geoms = idpp_interpolate(geoms, images_between=between)
        xyz_per_image = [geom.as_xyz() for geom in geoms]
        trj = "\n".join(xyz_per_image)
    # or just linear interpolation.
    elif between != 0:
        cos = ChainOfStates.ChainOfStates(geoms)
        cos.interpolate(between)
        geoms = cos.images
        xyz_per_image = [geom.as_xyz() for geom in geoms]
        trj = cos.as_xyz()

    return geoms


def dump_geoms(geoms, fn_base, trj_infix="", dump_trj=True):
    xyz_per_geom = [geom.as_xyz() for geom in geoms]
    if dump_trj:
        trj_str = "\n".join(xyz_per_geom)
        trj_fn = f"{fn_base}{trj_infix}.trj"
        with open(trj_fn, "w") as handle:
            handle.write(trj_str)
        print(f"Wrote all geometries to {trj_fn}.")
    for i, xyz in enumerate(xyz_per_geom):
        geom_fn = f"{fn_base}.geom_{i:03d}.xyz"
        with open(geom_fn, "w") as handle:
            handle.write(xyz)
        print(f"Wrote geom {i:03d} to {geom_fn}.")
    print()


def align(geoms):
    """Align all geometries onto the first using partical procrustes."""
    cos = ChainOfStates.ChainOfStates(geoms)
    procrustes(cos)
    return [geom for geom in cos.images]


def spline_redistribute(geoms):
    import numpy as np
    def get_coords_diffs(coords):
        cds = [0, ]
        for i in range(len(coords)-1):
            diff = np.linalg.norm(coords[i+1]-coords[i])
            cds.append(diff)
        cds = np.cumsum(cds)
        cds /= cds.max()
        return cds
    szts = SimpleZTS.SimpleZTS(geoms)
    pre_diffs = get_coords_diffs([image.coords for image in szts.images])
    szts.reparametrize()
    post_diffs = get_coords_diffs(szts.coords.reshape(-1,3))
    post_diffs = get_coords_diffs([image.coords for image in szts.images])
    cds_str = lambda cds: " ".join([f"{cd:.2f}" for cd in cds])
    print("Normalized path segments before splining:")
    print(cds_str(pre_diffs))
    print("Normalized path segments after redistribution along spline:")
    print(cds_str(post_diffs))
    return szts.images


def every(geoms, every_nth):
    # every_nth_geom = geoms[::every_nth]
    # The first geometry is always present, but the last geometry
    # may be missing.
    sampled_indices = list(range(0, len(geoms), every_nth))
    if sampled_indices[-1] != len(geoms)-1:
        sampled_indices.append(len(geoms)-1)
    sampled_inds_str = ", ".join([str(i) for i in sampled_indices])
    print(f"Sampled indices {sampled_inds_str}")
    # if every_nth_geom[-1] != geoms[-1]:
        # every_nth_geom.append(geoms[-1])
    every_nth_geom = [geoms[i] for i in sampled_indices]
    return every_nth_geom


def bohr2ang(xyzs):
    coords_angstrom = [geom.coords*0.529177249 for geom in geoms]
    raise Exception("Implement me")


def run():
    args = parse_args(sys.argv[1:])

    geoms = get_geoms(args.fns, comments=True)

    dump_trj = True
    trj_infix = ""
    if args.between:
        to_dump = get_geoms(args.fns, args.idpp, args.between)
        fn_base = "interpolated"
    elif args.align:
        to_dump = align(geoms)
        fn_base = "aligned"
    elif args.split:
        to_dump = geoms
        fn_base = "split"
        dump_trj = False
    elif args.reverse:
        to_dump = geoms[::-1]
        fn_base = "reversed"
    elif args.cleantrj:
        to_dump = geoms
        fn_base = "cleaned"
    elif args.first:
        to_dump = geoms[:args.first]
        fn_base = "first"
        trj_infix = f"_{args.first}"
    elif args.spline:
        to_dump = spline_redistribute(geoms)
        fn_base = "splined"
    elif args.every:
        to_dump = every(geoms, args.every)
        fn_base = "every"
        trj_infix = f"_{args.every}th"
    elif args.bohr2ang:
        to_dump = bohr2ang(args.xyz)
        fn_base = "angstrom"

    dump_geoms(to_dump, fn_base, trj_infix=trj_infix, dump_trj=dump_trj)

if __name__ == "__main__":
    run()
