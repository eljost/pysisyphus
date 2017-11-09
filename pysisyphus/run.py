#!/usr/bin/env python3

import argparse
import sys

from pysisyphus.calculators import *
from pysisyphus.cos import *
from pysisyphus.Geometry import Geometry
from pysisyphus.irc import *
from pysisyphus.optimizers import *
from pysisyphus.helpers import geom_from_xyz_file, geoms_from_trj

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

IRC_DICT = {
    "dvv": DampedVelocityVerlet.DampedVelocityVerlet,
    "euler": Euler.Euler,
    "gs": GonzalesSchlegel.GonzalesSchlegel,
    #"imk": IMKMod.IMKMod,
}


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("--cos", choices="neb szts".split(),
                        help="Chain of states method.")
    parser.add_argument("--opt", choices="bfgs fire sd".split(),
                        default="fire", help="Optimization algorithm.")
    parser.add_argument("--calc", choices="orca xtb".split(),
                        default="xtb", help="Calculator.")
    parser.add_argument("--addimgs", type=int,
                        default=8, help="Interpolate additional images.")
    parser.add_argument("--idpp", action="store_true")
    parser.add_argument("--align", action="store_true")
    parser.add_argument("--outdir")
    parser.add_argument("--irc", choices="".split())

    parser.add_argument("xyz", nargs="+")

    return parser.parse_args()


def run_cos(args):
    if len(args.xyz) == 1 and args.xyz[0].endswith(".trj"):
        geoms = geoms_from_trj(args.xyz[0])
    else:
        atoms_coords = [parse_xyz_file(fn) for fn in args.xyz]
        geoms = [Geometry(atoms, coords.flatten())
                 for atoms, coords in atoms_coords]
    """Handle this differently.
    Doing IDPP interpolation in a standalone calculator vs.
    doing linear interpolation in the COS object prevents
    a nice coherent treatment here. Maybe it's not that easy
    right now, because moving IDPP to ChainOfStates may lead
    to circular imports..."""
    if args.idpp:
        geoms = IDPP.idpp_interpolate(geoms, images_between=args.addimgs)
    cos = COS_DICT[args.cos](geoms)
    if not args.idpp and args.addimgs:
        cos.interpolate(args.addimgs)
    for image in cos.images:
        image.set_calculator(CALC_DICT[args.calc]())
    opt = OPT_DICT[args.opt](cos, align=args.align)
    opt.run()


def run_irc(args):
    geom = geom_from_xyz_file(args)
    irc = IRC_DICT[args.irc](geom)
    irc.run()
    #irc.write_trj(THIS_DIR, prefix)


def run_opt(args):
    geom = geom_from_xyz_file(args)
    opt = OPT_DICT[args.opt](geom)
    opt.run()


def run():
    args = parse_args(sys.argv[1:])

    # Do ChainOfStates method
    if args.cos:
        run_cos(args)
    # Do IRC
    elif args.irc:
        run_irc(args)
    # Do conventional optimization
    elif args.opt:
        run_opt(args)
    else:
        print("No run type specified!")

if __name__ == "__main__":
    run()
