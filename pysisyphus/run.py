#!/usr/bin/env python3

import argparse
import sys

from pysisyphus.calculators import *
from pysisyphus.cos import *
from pysisyphus.Geometry import Geometry
from pysisyphus.irc import *
from pysisyphus.optimizers import *
from pysisyphus.helpers import geom_from_xyz_file, geoms_from_trj


COS_DICT = {
    "neb": NEB.NEB,
    "szts": SimpleZTS.SimpleZTS,
}

CALC_DICT = {
    "orca": ORCA.ORCA,
    "xtb": XTB.XTB,
    "openmolcas": OpenMolcas.OpenMolcas
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

    parser.add_argument("--cos", choices=COS_DICT.keys(),
                        help="Chain of states method.")
    parser.add_argument("--opt", choices=OPT_DICT.keys(),
                        help="Optimization algorithm.")
    parser.add_argument("--calc", choices=CALC_DICT.keys(),
                        help="Calculator.")
    parser.add_argument("--interpol", type=int,
                        help="Interpolate additional images.")
    parser.add_argument("--idpp", action="store_true",
                        help="Use Image Dependent Pair Potential instead "
                        "of simple linear interpolation.")
    parser.add_argument("--align", action="store_true")
    parser.add_argument("--outdir")
    parser.add_argument("--irc", choices=IRC_DICT.keys())
    parser.add_argument("--parallel", type=int, default=0)

    parser.add_argument("xyz", nargs="+")

    return parser.parse_args()


def get_geoms(args):
    # Read .xyz or .trj files
    if len(args.xyz) == 1 and args.xyz[0].endswith(".trj"):
        geoms = geoms_from_trj(args.xyz[0])
    else:
        geoms = [geom_from_xyz_file(fn) for fn in args.xyz]

    # Do IDPP interpolation if requested,
    if args.idpp:
        geoms = IDPP.idpp_interpolate(geoms, images_between=args.interpol)
    # or just linear interpolation.
    elif args.interpol:
        cos = ChainOfStates.ChainOfStates(geoms)
        cos.interpolate(args.interpol)
        geoms = cos.images

    return geoms


def run_cos(args):
    geoms = get_geoms(args)
    cos = COS_DICT[args.cos](geoms, parallel=args.parallel)
    for i, image in enumerate(cos.images):
        name = f"image_{i:03d}"
        image.set_calculator(CALC_DICT[args.calc](name=name))
    kwargs = {
        #"max_cycles": 2,
    }
    opt = OPT_DICT[args.opt](cos, align=args.align, **kwargs)
    opt.run()


def run_irc(args):
    assert(len(arg.xyz) == 1)
    geom = get_geoms(args)[0]
    geom.set_calculator(CALC_DICT[args.calc]())
    irc = IRC_DICT[args.irc](geom)
    irc.run()
    #irc.write_trj(THIS_DIR, prefix)


def run_opt(args):
    assert(len(args.xyz) == 1)
    geom = get_geoms(args)[0]
    geom.set_calculator(CALC_DICT[args.calc]())
    opt = OPT_DICT[args.opt](geom)
    opt.run()


def run_interpolation(args):
    geoms = get_geoms(args)
    trj_fn = "interpolated.trj"
    trj_str = "\n".join([geom.as_xyz() for geom in geoms])
    with open(trj_fn, "w") as handle:
        handle.write(trj_str)
    print(f"Wrote interpolated coordinates to {trj_fn}.")


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
    # Just interpolate
    elif args.interpol:
        run_interpolation(args)
    else:
        print("Please specify a run type!")

if __name__ == "__main__":
    run()
