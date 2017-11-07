#!/usr/bin/env python3

import argparse
import sys

from pysisyphus.calculators import *
from pysisyphus.cos import *
from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers import *

from qchelper.geometry import parse_xyz_file

COS_DICT = {
    "neb": NEB.NEB,
    "szts": SimpleZTS.SimpleZTS,
}

CALC_DICT = {
    "orca": ORCA.ORCA,
    "xtb": XTB.XTB,
}

OPT_DICT = {
    "fire": FIRE.FIRE,
    "bfgs": BFGS.BFGS,
    "sd": SteepestDescent.SteepestDescent,
}

def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("--cos", choices="neb szts".split(),
                        default="neb")
    parser.add_argument("--opt", choices="bfgs fire sd".split(),
                        default="fire")
    parser.add_argument("--calc", choices="orca xtb".split(),
                        default="xtb")
    parser.add_argument("--addimgs", type=int,
                        default=8)
    parser.add_argument("--idpp", action="store_true")
    parser.add_argument("--outdir")

    parser.add_argument("xyz", nargs="+")

    return parser.parse_args()


def run_cos(args):
    atoms_coords = [parse_xyz_file(fn) for fn in args.xyz]
    geoms = [Geometry(atoms, coords.flatten())
             for atoms, coords in atoms_coords]
    """Handle this differently.
    Doing IDPP interpolation in a standalone calculator vs.
    doing linear interpolation in the COS object prevents
    an easy treatment here."""
    if args.idpp:
        geoms = IDPP.idpp_interpolate(geoms, images_between=args.addimgs)
    cos = COS_DICT[args.cos](geoms)
    if not args.idpp:
        cos.interpolate(args.addimgs)
    for image in cos.images:
        image.set_calculator(CALC_DICT[args.calc]())
    #kwargs["out_dir"] = "neb_fire"
    opt = OPT_DICT[args.opt](cos)#, **kwargs)
    opt.run()


def run():
    args = parse_args(sys.argv[1:])

    if args.cos:
        run_cos(args)

if __name__ == "__main__":
    run()
