#!/usr/bin/env python3

import argparse
import sys

import yaml

from pysisyphus.calculators import *
from pysisyphus.cos import *
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_xyz_file, geoms_from_trj
from pysisyphus.irc import *
from pysisyphus.optimizers import *


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
    parser.add_argument("--xyz", nargs="+")
    parser.add_argument("--yaml")

    return parser.parse_args()


def get_geoms(xyz_fns, idpp=False, between=0):
    # Read .xyz or .trj files
    if len(xyz_fns) == 1 and args.xyz[0].endswith(".trj"):
        geoms = geoms_from_trj(args.xyz[0])
    else:
        geoms = [geom_from_xyz_file(fn) for fn in xyz_fns]

    # Do IDPP interpolation if requested,
    if idpp:
        geoms = IDPP.idpp_interpolate(geoms, images_between=between)
    # or just linear interpolation.
    elif interpol:
        cos = ChainOfStates.ChainOfStates(geoms)
        cos.interpolate(between)
        geoms = cos.images

    return geoms


"""
def run_cos(args):
    geoms = get_geoms(args)
    cos = COS_DICT[args.cos](geoms, parallel=args.parallel)
    for i, image in enumerate(cos.images):
        name = f"image_{i:03d}"
        image.set_calculator(CALC_DICT[args.calc](name=name))
    #opt = OPT_DICT[args.opt](cos, align=args.align, **kwargs)
    opt = OPT_DICT[args.opt](cos, align=args.align)
    opt.run()
"""

def run_cos(cos, get_calc, get_opt):
    for i, image in enumerate(cos.images):
        name = f"image_{i:03d}"
        image.set_calculator(get_calc({"name":name}))
    opt = get_opt(cos)
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


def get_defaults(conf_dict):
    # dd = default_dict
    dd = dict()
    if "cos" in conf_dict:
        dd["cos"] = {
            "type": "neb",
            "parallel": 0,
        }
        dd["opt"] = {
            "type": "bfgs",
            "align": True,
        }
        dd["interpol"] = {
            "idpp": True,
            "between": 8,
        }

    return dd


def handle_yaml(yaml_str):
    yaml_dict = yaml.load(yaml_str)
    run_dict = get_defaults(yaml_dict)
    # Update nested entries
    key_set = set(yaml_dict.keys())
    for key in key_set & set(("cos", "opt", "interpol")):
        run_dict[key].update(yaml_dict[key])
    # Update non nested entries
    for key in key_set & set(("calc", "charge", "mult", "xyz")):
        run_dict[key] = yaml_dict[key]
    print(run_dict)

    xyz = run_dict["xyz"]
    if run_dict["interpol"]:
        idpp = run_dict["interpol"]["idpp"]
        between = run_dict["interpol"]["between"]
    if run_dict["opt"]:
        opt_key = run_dict["opt"].pop("type")
        opt_kwargs = run_dict["opt"]
    if run_dict["cos"]:
        cos_key = run_dict["cos"].pop("type")
        cos_kwargs = run_dict["cos"]

    calc_key = run_dict["calc"].pop("type")
    calc_kwargs = run_dict["calc"]
    get_calc = lambda kw: CALC_DICT[calc_key](**kw, **calc_kwargs)
    get_opt = lambda geoms: OPT_DICT[opt_key](geoms, **opt_kwargs)

    geoms = get_geoms(xyz, idpp, between)
    if run_dict["cos"]:
        cos = COS_DICT[cos_key](geoms, **cos_kwargs)
        run_cos(cos, get_calc, get_opt)



def run():
    args = parse_args(sys.argv[1:])

    # Do ChainOfStates method
    if args.yaml:
        with open(args.yaml) as handle:
            yaml_str = handle.read()
        handle_yaml(yaml_str)
    elif args.cos:
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
