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


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("--idpp", action="store_true",
                        help="Use Image Dependent Pair Potential instead "
                             "of simple linear interpolation.")
    parser.add_argument("--xyz", nargs="+")

    action_group = parser.add_mutually_exclusive_group(required=True)
    action_group.add_argument("--between", type=int,
                              help="Interpolate additional images.")
    action_group.add_argument("--align", nargs="+",
                              help="Align geometries onto the first geometry "
                                   "read from multiple .xyz or one .trj file.")
    action_group.add_argument("--split",
                              help="Split a supplied .trj file in multiple "
                                   ".xyz files.")
    action_group.add_argument("--reverse",
                              help="Reverse a .trj file.")
    action_group.add_argument("--cleantrj",
                              help="Clean a .trj file.")
    action_group.add_argument("--spline",
                              help="Evenly redistribute geometries along a "
                                   "splined path.")
    return parser.parse_args()


def get_geoms(xyz_fns, idpp=False, between=0, dump=False, multiple_geoms=False):
    """Returns a list of Geometry objects."""
    # Handle a single .xyz file
    if isinstance(xyz_fns, str) and xyz_fns.endswith(".xyz"):
        geoms = [geom_from_xyz_file(xyz_fns), ]
    # Handle a single .trj file
    elif len(xyz_fns) == 1 and xyz_fns[0].endswith(".trj"):
        geoms = geoms_from_trj(xyz_fns[0])
    # How is this different from above?
    elif isinstance(xyz_fns, str) and xyz_fns.endswith(".trj"):
        geoms = geoms_from_trj(xyz_fns)
    elif multiple_geoms:
        geoms = geoms_from_trj(xyz_fns[0])
    # Handle multiple .xyz files
    else:
        geoms = [geom_from_xyz_file(fn) for fn in xyz_fns]

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

    if dump and len(geoms) > 2:
        dump_geometry_strings("interpolated", trj, xyz_per_image)

    return geoms


def dump_geometry_strings(base, trj="", xyz_per_image=[]):
    if trj:
        trj_fn = f"{base}.trj"
        with open(trj_fn, "w") as handle:
            handle.write(trj)
        print(f"Wrote all geometries to {trj_fn}.")
    for i, xyz in enumerate(xyz_per_image):
        image_fn = f"{base}.image_{i}.xyz"
        with open(image_fn, "w") as handle:
            handle.write(xyz)
        print(f"Wrote image {i} to {image_fn}.")
    print()


def run_interpolation(args):
    """Interpolate between given geometries."""
    geoms = get_geoms(args.xyz, args.idpp, args.between, dump=True)


def align(fns):
    """Align all geometries onto the first using partical procrustes."""
    geoms = get_geoms(fns, multiple_geoms=True)
    cos = ChainOfStates.ChainOfStates(geoms)
    procrustes(cos)
    trj = cos.as_xyz()
    xyz_per_image = [image.as_xyz() for image in cos.images]
    dump_geometry_strings("aligned", trj, xyz_per_image)


def split(trj_fn):
    """Split a .trj in several .xyz files."""
    geoms = get_geoms(trj_fn, multiple_geoms=True)
    xyz_per_image = [geom.as_xyz() for geom in geoms]
    dump_geometry_strings("split", xyz_per_image=xyz_per_image)


def reverse_trj(trj_fn):
    """Reverse the ordering of the geometries in a .trj file."""
    geoms = get_geoms(trj_fn, multiple_geoms=True)
    xyz_per_image = list(reversed([geom.as_xyz() for geom in geoms]))
    reversed_trj = "\n".join(xyz_per_image)
    dump_geometry_strings("reversed", trj=reversed_trj)


def clean_trj(trj_fn):
    """Drops everything from an extended .trj file, e.g. additional columns
    containing gradients and so on. Keeps only the first four columns
    specifying the atom and xyz coordinates."""
    geoms = get_geoms(trj_fn, multiple_geoms=True)
    xyz_per_image = [geom.as_xyz() for geom in geoms]
    trj = "\n".join(xyz_per_image)
    dump_geometry_strings("cleaned", trj=trj)


def spline(trj_fn):
    import numpy as np
    def get_coords_diffs(coords):
        cds = [0, ]
        for i in range(len(coords)-1):
            diff = np.linalg.norm(coords[i+1]-coords[i])
            cds.append(diff)
        cds = np.cumsum(cds)
        cds /= cds.max()
        return cds
    geoms = get_geoms(trj_fn, multiple_geoms=True)
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
    xyz_per_image = [image.as_xyz() for image in szts.images]
    trj = "\n".join(xyz_per_image)
    dump_geometry_strings("splined", trj=trj)


def run():
    args = parse_args(sys.argv[1:])

    if args.between:
        run_interpolation(args)
    elif args.align:
        align(args.align)
    elif args.split:
        split(args.split)
    elif args.reverse:
        reverse_trj(args.reverse)
    elif args.cleantrj:
        clean_trj(args.cleantrj)
    elif args.spline:
        spline(args.spline)

if __name__ == "__main__":
    run()
