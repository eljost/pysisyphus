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

from pysisyphus.calculators import *
from pysisyphus.cos import *
from pysisyphus.irc import *
from pysisyphus.optimizers import *
from pysisyphus.trj import get_geoms


COS_DICT = {
    "neb": NEB.NEB,
    "szts": SimpleZTS.SimpleZTS,
}

CALC_DICT = {
    "orca": ORCA.ORCA,
    "xtb": XTB.XTB,
    "openmolcas": OpenMolcas.OpenMolcas,
    "g09": Gaussian09.Gaussian09,
    "turbomole": Turbomole.Turbomole,
}

OPT_DICT = {
    "fire": FIRE.FIRE,
    # Removing BFGS for now until save_also is implemented
    # and rotating the hessian works properly
    "bfgs": BFGS.BFGS,
    "sd": SteepestDescent.SteepestDescent,
    "cg": ConjugateGradient.ConjugateGradient,
    "qm": QuickMin.QuickMin,
    "scipy": SciPyOptimizer.SciPyOptimizer,
    "rfo": RFOptimizer.RFOptimizer,
}

IRC_DICT = {
    "dvv": DampedVelocityVerlet.DampedVelocityVerlet,
    "euler": Euler.Euler,
    "gs": GonzalesSchlegel.GonzalesSchlegel,
    #"imk": IMKMod.IMKMod,
}


def parse_args(args):
    parser = argparse.ArgumentParser()

    action_group = parser.add_mutually_exclusive_group(required=True)
    action_group.add_argument("yaml", nargs="?",
                              help="Start pysisyphus with input from a "
                                   "YAML file.")
    action_group.add_argument("--clean", action="store_true")

    parser.add_argument("--restart", action="store_true",
                        help="Continue a previously crashed/aborted/... "
                             "pysisphus run.")
    return parser.parse_args()


def get_calc(index, base_name, calc_key, calc_kwargs):
    # Expand values that contain the $IMAGE pattern over all images.
    # We have to use a copy of calc_kwargs to keep the $IMAGE pattern.
    # Otherwise it would be replace at it's first occurence and would
    # be gone in the following items.
    kwargs_copy = copy.deepcopy(calc_kwargs)
    for key, value in kwargs_copy.items():
        if not isinstance(value, str) or not ("$IMAGE" in value):
            continue
        kwargs_copy[key] = value.replace("$IMAGE", f"{index:03d}")
    kwargs_copy["base_name"] = base_name
    kwargs_copy["calc_number"] = index
    return CALC_DICT[calc_key](**kwargs_copy)



def run_cos(cos, calc_getter, opt_getter):
    for i, image in enumerate(cos.images):
        image.set_calculator(calc_getter(i))
    opt = opt_getter(cos)
    opt.run()


def run_opt(geom, calc_getter, opt_getter):
    geom.set_calculator(calc_getter(0))
    opt = opt_getter(geom)
    opt.run()


"""
def run_irc(args):
    assert(len(arg.xyz) == 1)
    geom = get_geoms(args)[0]
    geom.set_calculator(CALC_DICT[args.calc]())
    irc = IRC_DICT[args.irc](geom)
    irc.run()
    #irc.write_trj(THIS_DIR, prefix)
"""


def get_defaults(conf_dict):
    # Defaults
    dd = {
        "interpol": {
            "idpp": False,
            "between": 0,
        },
        "cos": None,
    }
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
    elif "opt" in conf_dict:
        dd["opt"] = {
            "type": "cg",
            "dump": True,
            "alpha": 0.25,
        }

    return dd


def get_last_calc_cycle():
    def keyfunc(path):
        return re.match("image_\d+.(\d+).out", str(path))[1]
    cwd = Path(".")
    calc_logs = [str(cl) for cl in cwd.glob("image_*.*.out")]
    calc_logs = sorted(calc_logs, key=keyfunc)
    grouped = itertools.groupby(calc_logs, key=keyfunc)
    # Find the last completly finished cycle.
    last_length = 0
    last_calc_cycle = 0
    for calc_cycle, group in grouped:
        cycle_length = len(list(group))
        if cycle_length < last_length:
            # When this is True we have a cycle that has less
            # items than last one, that is an unfinished cycle.
            break
        last_length = cycle_length
        last_calc_cycle = int(calc_cycle)
    if last_calc_cycle == 0:
        print("Can't find any old calculator logs.")
    print(f"Last calculation counter is {last_calc_cycle}.")
    return last_calc_cycle


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
    return run_dict


def main(run_dict, restart):
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

    if restart:
        print("Trying to restart calculation. Skipping interpolation.")
        idpp = False
        between = 0
        # Load geometries of latest cycle
        cwd = Path(".")
        trjs = [str(trj) for trj in cwd.glob("cycle_*.trj")]
        if len(trjs) == 0:
            print("Can't restart. Found no previous coordinates.")
            sys.exit()
        xyz = natsorted(trjs)[-1]
        last_cycle = int(re.search("(\d+)", xyz)[0])
        print(f"Last cycle was {last_cycle}.")
        print(f"Using '{xyz}' as input geometries.")
        opt_kwargs["last_cycle"] = last_cycle
        last_calc_cycle = get_last_calc_cycle()
        run_dict["calc"]["last_calc_cycle"] = last_calc_cycle

    calc_key = run_dict["calc"].pop("type")
    calc_kwargs = run_dict["calc"]
    calc_getter = lambda index: get_calc(index, "image", calc_key, calc_kwargs)
    opt_getter = lambda geoms: OPT_DICT[opt_key](geoms, **opt_kwargs)

    geoms = get_geoms(xyz, idpp, between, dump=True)
    if run_dict["cos"]:
        cos = COS_DICT[cos_key](geoms, **cos_kwargs)
        run_cos(cos, calc_getter, opt_getter)
    elif run_dict["opt"]:
        assert(len(geoms) == 1)
        run_opt(geoms[0], calc_getter, opt_getter)


def clean():
    """Deletes files from previous runs in the cwd.
    A similar function could be used to store everything ..."""
    cwd = Path(".").resolve()
    rm_globs = (
        "image*.trj",
        "image*.out",
        "cycle*.trj",
        "interpolated.trj",
        "interpolated.image*.xyz",
        "calculator.log",
        "optimizer.log",
        "optimization.trj",
        "cos.log",
        "*.gradient",
        "image_results.yaml",
        "optimizer_results.yaml",
        # ORCA specific
        "image*.gbw",
        "image*.engrad",
        "image*.hessian",
        # OpenMOLCAS specific
        "image*.RasOrb",
        "image*.in",
        "image*.JobIph",
        "calculator*.out",
        "*rasscf.molden",
        # Gaussian specific
        "image*.fchk",
        "image*.log",
        # Turbomole specific
        "image*.mos",
        "image*.alpha",
        "image*.beta",
        "image*.control",
        "image*.ciss_a",
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

    if args.yaml:
        with open(args.yaml) as handle:
            yaml_str = handle.read()
        run_dict = handle_yaml(yaml_str)
        pprint(run_dict)
        main(run_dict, args.restart)
    elif args.clean:
        clean()

if __name__ == "__main__":
    run()
