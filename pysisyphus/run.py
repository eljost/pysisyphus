#!/usr/bin/env python3

import argparse
import os
from pathlib import Path
from pprint import pprint
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
    "cg": ConjugateGradient.ConjugateGradient,
    "qm": QuickMin.QuickMin,
    "scipy": SciPyOptimizer.SciPyOptimizer,
}

IRC_DICT = {
    "dvv": DampedVelocityVerlet.DampedVelocityVerlet,
    "euler": Euler.Euler,
    "gs": GonzalesSchlegel.GonzalesSchlegel,
    #"imk": IMKMod.IMKMod,
}


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("--between", type=int,
                        help="Interpolate additional images.")
    parser.add_argument("--idpp", action="store_true",
                        help="Use Image Dependent Pair Potential instead "
                        "of simple linear interpolation.")
    parser.add_argument("--xyz", nargs="+")
    parser.add_argument("--yaml")
    parser.add_argument("--clean", action="store_true")

    return parser.parse_args()


def get_calc(index, name_base, calc_key, calc_kwargs):
    """
    import copy
    kwargs_copy = copy.copy(calc_kwargs)
    for key in kwargs_copy:
        val = kwargs_copy[key]
        if "$IMAGE" in str(val):
            kwargs_copy[key] = val.replace("$IMAGE", str(index))
    kwargs_copy["name"] = f"{name_base}_{index:03d}"
    """
    kwargs = {key: str(calc_kwargs[key]).replace("$IMAGE", "{index:03d}")
              for key in calc_kwargs}
    kwargs["name"] = f"{name_base}_{index:03d}"
    return CALC_DICT[calc_key](**kwargs)


def get_geoms(xyz_fns, idpp=False, between=0):
    # Read .xyz or .trj files
    if len(xyz_fns) == 1 and args.xyz[0].endswith(".trj"):
        geoms = geoms_from_trj(args.xyz[0])
    else:
        geoms = [geom_from_xyz_file(fn) for fn in xyz_fns]

    # Do IDPP interpolation if requested,
    trj = ""
    xyz_per_image = list()
    if idpp:
        geoms = IDPP.idpp_interpolate(geoms, images_between=between)
        xyz_per_image = [geom.as_xyz() for geom in geoms]
        trj = "\n".join(xyz_per_image)
    # or just linear interpolation.
    elif between != 0:
        cos = ChainOfStates.ChainOfStates(geoms)
        cos.interpolate(between)
        geoms = cos.images
        xyz_per_image = [geom.as_xyz() for geom in geoms]
        trj = cos.as_xyz()

    if trj:
        trj_fn = "interpolated.trj"
        with open(trj_fn, "w") as handle:
            handle.write(trj)
        print(f"Wrote interpolated coordinates to {trj_fn}.")
    for i, xyz in enumerate(xyz_per_image):
        image_fn = f"interpolated.image_{i}.xyz"
        with open(image_fn, "w") as handle:
            handle.write(xyz)
        print(f"Wrote image {i} to {image_fn}.")

    return geoms


def run_cos(cos, calc_getter, get_opt):
    for i, image in enumerate(cos.images):
        image.set_calculator(calc_getter(i))
    opt = get_opt(cos)
    opt.run()


"""
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
"""


def run_interpolation(args):
    geoms = get_geoms(args.xyz, args.idpp, args.between)
    """
    trj_fn = "interpolated.trj"
    trj_str = "\n".join([geom.as_xyz() for geom in geoms])
    with open(trj_fn, "w") as handle:
        handle.write(trj_str)
    """


def get_defaults(conf_dict):
    # dd = default_dict
    dd = dict()
    if "cos" in conf_dict:
        dd["cos"] = {
            "type": "neb",
            "parallel": 0,
        }
        dd["opt"] = {
            "type": "cg",
            "align": True,
            "dump": True,
        }
        dd["interpol"] = {
            "idpp": False,
            "between": 0,
        }

    return dd


def handle_yaml(yaml_str):
    yaml_dict = yaml.load(yaml_str)
    # Load defaults to have a sane baseline
    run_dict = get_defaults(yaml_dict)
    # Update nested entries
    key_set = set(yaml_dict.keys())
    for key in key_set & set(("cos", "opt", "interpol")):
        run_dict[key].update(yaml_dict[key])
    # Update non nested entries
    for key in key_set & set(("calc", "xyz")):
        run_dict[key] = yaml_dict[key]
    pprint(run_dict)

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
    calc_getter = lambda index: get_calc(index, "image", calc_key, calc_kwargs)
    get_opt = lambda geoms: OPT_DICT[opt_key](geoms, **opt_kwargs)

    geoms = get_geoms(xyz, idpp, between)
    if run_dict["cos"]:
        cos = COS_DICT[cos_key](geoms, **cos_kwargs)
        run_cos(cos, calc_getter, get_opt)


def clean():
    """Deletes files from previous runs in the cwd.
    A similar function could be used to store everything ..."""
    cwd = Path(".").resolve()
    rm_globs = (
        "image*.trj",
        "image*.out",
        "cycle*.trj",
        "results.yaml",
        "interpolated.trj",
        "interpolated.image*.xyz",
        "calculator.log",
        "optimizer.log",
        "optimized.trj",
        "cos.log",
        # ORCA specific
        "image*.gbw",
        "image*.engrad",
        "image*.hessian",
        # OpenMOLCAS specific
        "image*.RasOrb",
    )
    to_rm_paths = list()
    for glob in rm_globs:
        to_rm_paths.extend(list(cwd.glob(glob)))
    to_rm_strs = [str(p) for p in to_rm_paths]
    for s in to_rm_strs:
        print(s)
    rly_delete = input("Delete these files? (yes/no)\n")
    if rly_delete != "yes":
        print("Aborting")
        return
    else:
        for p in to_rm_paths:
            os.remove(p)
            print(f"Deleted {p}")

def run():
    args = parse_args(sys.argv[1:])

    # Do ChainOfStates method
    if args.yaml:
        with open(args.yaml) as handle:
            yaml_str = handle.read()
        handle_yaml(yaml_str)
    elif args.between:
        run_interpolation(args)
    elif args.clean:
        clean()
    else:
        print("Please specify a run type! Show help with -h.")

if __name__ == "__main__":
    run()
