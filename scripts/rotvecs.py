#!/usr/bin/env python3

import argparse
import sys

import numpy as np

from pysisyphus.helpers import geom_loader
from pysisyphus.trj import align_coords3d_onto_vec


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("fn")
    parser.add_argument("inds", nargs=2, type=int)
    parser.add_argument("ref_vec", nargs=3, type=float)
    parser.add_argument("--jmol", action="store_true")

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    fn = args.fn
    inds = args.inds
    ref_vec = np.array(args.ref_vec)
    geom = geom_loader(fn)

    # Rotates coordiantes and new geometry to store them
    c3dr = align_coords3d_onto_vec(geom.coords3d, *inds, ref_vec)
    geom_rot = geom.copy()
    geom_rot.coords3d = c3dr
    xyz_rot = geom.as_xyz()
    print(f"Rotated coordinates:\n{xyz_rot}")

    if args.jmol:
        geom_rot.jmol()


if __name__ == "__main__":
    run()
