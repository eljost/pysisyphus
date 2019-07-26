#!/usr/bin/env python3

import argparse
import cloudpickle
import copy
import itertools
import os
from pathlib import Path
from pprint import pprint
import re
import shutil
import sys
import time

from distributed import Client
from natsort import natsorted
import numpy as np
import yaml

from pysisyphus.calculators import *
from pysisyphus.cos import *
from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.cos.GrowingChainOfStates import GrowingChainOfStates
from pysisyphus.overlaps.Overlapper import Overlapper
from pysisyphus.overlaps.couplings import couplings
from pysisyphus.overlaps.sorter import sort_by_overlaps
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_xyz_file, confirm_input, shake_coords
from pysisyphus.irc import *
from pysisyphus.stocastic import *
from pysisyphus.init_logging import init_logging
from pysisyphus.optimizers import *
from pysisyphus.tsoptimizers import *
from pysisyphus.trj import get_geoms, dump_geoms
from pysisyphus.tsoptimizers.dimer import dimer_method
from pysisyphus._version import get_versions


CALC_DICT = {
    "g09": Gaussian09.Gaussian09,
    "g16": Gaussian16.Gaussian16,
    "openmolcas": OpenMolcas.OpenMolcas,
    "orca": ORCA.ORCA,
    "psi4": Psi4,
    # "pyscf": PySCF,
    "turbomole": Turbomole.Turbomole,
    "xtb": XTB.XTB,
}

try:
    from pysisyphus.calculators.PySCF import PySCF
    CALC_DICT["pyscf"] = PySCF
except ModuleNotFoundError:
    pass

COS_DICT = {
    "neb": NEB.NEB,
    "aneb": AdaptiveNEB.AdaptiveNEB,
    "feneb": FreeEndNEB.FreeEndNEB,
    "szts": SimpleZTS.SimpleZTS,
    "gs": GrowingString.GrowingString,
    "fs": FreezingString.FreezingString,
}

from pysisyphus.optimizers.BFGS_ import BFGS_
OPT_DICT = {
    "fire": FIRE.FIRE,
    # Removing BFGS for now until save_also is implemented
    # and rotating the hessian works properly
    "bfgs": BFGS.BFGS,
    "lbfgs": LBFGS.LBFGS,
    #"lbfgsm": LBFGS_mod.LBFGS,
    "sd": SteepestDescent.SteepestDescent,
    "cg": ConjugateGradient.ConjugateGradient,
    "qm": QuickMin.QuickMin,
    "scipy": SciPyOptimizer.SciPyOptimizer,
    "rfo": RFOptimizer.RFOptimizer,
    "rsrfo": RSRFOptimizer,
    "rsa": RSAlgorithm.RSAlgorithm,
    "anc": ANCOptimizer.ANCOptimizer,
    "string": StringOptimizer.StringOptimizer,
    "sqnm": StabilizedQNMethod.StabilizedQNMethod,
    "bfgs_": BFGS_,
}

TSOPT_DICT = {
    "dimer": dimer_method,
    "rsprfo": RSPRFOptimizer,
}

IRC_DICT = {
    # "dvv": DampedVelocityVerlet.DampedVelocityVerlet,
    "euler": Euler.Euler,
    "gs": GonzalesSchlegel.GonzalesSchlegel,
    "imk": IMKMod.IMKMod,
    "rk4": RK4,
}

STOCASTIC_DICT = {
    "frag": FragmentKick.FragmentKick,
    "kick": Kick.Kick,
}


def parse_args(args):
    parser = argparse.ArgumentParser()

    action_group = parser.add_mutually_exclusive_group(required=True)
    action_group.add_argument("yaml", nargs="?",
        help="Start pysisyphus with input from a YAML file."
    )
    action_group.add_argument("--clean", action="store_true",
        help="Ask for confirmation before cleaning."
    )
    action_group.add_argument("--fclean", action="store_true",
        help="Force cleaning without prior confirmation."
    )
    action_group.add_argument("--couplings", type=int, nargs="+",
        help="Create coupling elements."
    )
    action_group.add_argument("--sort-by-overlaps", nargs=2,
                              metavar=("ENERGIES", "OVERLAPS"),
        help="Sort precomputed energies of (excited) states by precomputed "
             "overlaps."
    )

    run_type_group = parser.add_mutually_exclusive_group(required=False)
    run_type_group.add_argument("--dryrun", action="store_true",
        help="Only generate a sample input (if meaningful) for checking."
    )
    run_type_group.add_argument("--restart", action="store_true",
        help="Continue a previously crashed/aborted/... pysisphus run."
    )
    run_type_group.add_argument("--overlaps", action="store_true",
        help="Calculate overlaps between transition density matrices "
             "(tden) or wavefunctions (wf)."
    )
    run_type_group.add_argument("--cp", "--copy",
        help="Copy .yaml file and corresponding geometries to a "
             "new directory. Similar to TURBOMOLEs cpc command."
    )

    parser.add_argument("--scheduler", default=None,
        help="Address of the dask scheduler."
    )
    parser.add_argument("--consider-first", type=int, default=None,
        help="Consider the first N states for sorting according "
             "to overlaps."
    )

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


def get_calc_closure(base_name, calc_key, calc_kwargs):
    index = 0
    def calc_getter():
        nonlocal index
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
        index += 1
        return CALC_DICT[calc_key](**kwargs_copy)
    return calc_getter


def preopt_ends(xyz, calc_getter):
    """Run optimization on first and last geometry in xyz and return
    updated xyz variable containing the optimized ends and any
    intermediate image that was present in the original list."""
    geoms = get_geoms(xyz, coord_type="redund")
    assert len(geoms) >= 2, "Need at least two geometries!"

    # middle = geoms[1:-1]

    def opt_getter(geom):
        opt_kwargs = {
            "max_cycles": 150,
            "thresh": "gau",
            "trust_max": 0.3,
        }
        opt = RFOptimizer.RFOptimizer(geom, **opt_kwargs)
        return opt

    out_xyz = list()
    for ind, str_ in ((0, "first"), (-1, "last")):
        print(f"Preoptimizing {str_} geometry.")
        geom = geoms[ind]
        opt = run_opt(geom, calc_getter, opt_getter)
        if not opt.is_converged:
            print("Problem in preoptimization of {str_}. Exiting!")
            sys.exit()
        print(f"Preoptimization of {str_} geometry converged!")
        opt_fn = f"{str_}_preopt.xyz"
        shutil.move("final_geometry.xyz", opt_fn)
        print(f"Saved final preoptimized structure to '{opt_fn}'.")
        out_xyz.append(opt_fn)
        print()
    return out_xyz


def run_cos(cos, calc_getter, opt_getter):
    for i, image in enumerate(cos.images):
        image.set_calculator(calc_getter(i))
    opt = opt_getter(cos)
    opt.run()
    hei_coords, hei_energy, hei_tangent = cos.get_splined_hei()
    hei_geom = Geometry(cos.images[0].atoms, hei_coords)
    hei_geom.coords = hei_coords
    hei_fn = "splined_hei.xyz"
    with open(hei_fn, "w") as handle:
        handle.write(hei_geom.as_xyz())
    print(f"Wrote splined HEI to '{hei_fn}'")


def run_cos_tsopt(cos, tsopt_key, tsopt_kwargs, calc_getter=None):
    print("Starting TS optimization after chain of states method.")

    # Use plain HEI
    hei_index = cos.get_hei_index()
    hei_image = cos.images[hei_index]
    prim_indices = (None if cos.coord_type == "cart"
                    else hei_image.internal.prim_indices
    )
    ts_geom = Geometry(hei_image.atoms, hei_image.cart_coords,
                       coord_type="redund", prim_indices=prim_indices)

    hei_tangent = cos.get_tangent(hei_index)
    # Convert tangent to redundant internals
    if cos.coord_type == "dlc":
        redund_tangent = hei_image.internal.Ut_inv @ hei_tangent
    elif cos.coord_type == "cart":
        redund_tangent = ts_geom.internal.B_prim @ hei_tangent
    else:
        raise Exception("Unknown coord_type!")

    # Use splined HEI
    # hei_coords, hei_energy, hei_tangent = cos.get_splined_hei()
    # ts_geom = Geometry(cos.images[0].atoms, hei_coords, coord_type="redund")
    ts_geom.set_calculator(calc_getter())
    print("Splined HEI (TS guess)")
    print(ts_geom.as_xyz())
    hei_xyz_fn = "splined_hei.xyz"
    with open(hei_xyz_fn, "w") as handle:
        handle.write(ts_geom.as_xyz())
    print(f"Wrote splined HEI coordinates to '{hei_xyz_fn}'")
    hei_tangent_fn = "splined_hei_tangent"
    np.savetxt("hei_tangent", hei_tangent)
    print(f"Wrote splined HEI tangent to '{hei_tangent_fn}'")

    ts_optimizer = TSOPT_DICT[tsopt_key]
    do_hess = tsopt_kwargs.pop("do_hess")

    if tsopt_key == "dimer":
        geoms = [ts_geom, ]
        tsopt_kwargs.update({
            "restrict_step": "max",
            "N_init": hei_tangent,
            "calc_getter": calc_getter,
        })
        dimer_result = ts_optimizer(geoms, **tsopt_kwargs)
        dimer_cycles = dimer_result.dimer_cycles
        last_cycle = dimer_cycles[-1]
        ts_coords = last_cycle.trans_coords[1]
    elif tsopt_key == "rsprfo":
        # # Determine which imaginary mode has the highest overlap
        # # with the splined HEI tangent.
        eigvals, eigvecs = np.linalg.eigh(ts_geom.hessian)
        neg_inds = eigvals < 0
        eigval_str = np.array2string(eigvals[neg_inds], precision=6)
        print(f"Negative eigenvalues at splined HEI: {eigval_str}")
        neg_eigvecs = eigvecs.T[neg_inds]
        ovlps = [np.abs(imag_mode.dot(redund_tangent)) for imag_mode in neg_eigvecs]
        print("Overlaps between HEI tangent and imaginary modes:")
        for i, ov in enumerate(ovlps):
            print(f"\t{i:02d}: {ov:.6f}")
        root = np.argmax(ovlps)
        print(f"Imaginary mode {root} has highest overlap with splined HEI tangent")
        # Use mode with highest overlap as initial root
        tsopt_kwargs["root"] = root
        prfo = ts_optimizer(ts_geom, **tsopt_kwargs)
        prfo.run()

    print(f"Optimized TS coords:")
    print(ts_geom.as_xyz())
    with open("ts_opt.xyz", "w") as handle:
        handle.write(ts_geom.as_xyz())

    if do_hess:
        hessian = ts_geom.hessian
        eigvals, eigvecs = np.linalg.eigh(hessian)
        ev_thresh = -1e-4
        neg_eigvals = eigvals[eigvals < ev_thresh]
        print(f"Self found {neg_eigvals.size} eigenvalues < {ev_thresh}.")
        print("Negative eigenvalues: ", neg_eigvals)


# def run_cos_dimer(cos, dimer_kwargs, calc_getter):
    # print("Starting dimer method")
    # hei_coords, hei_energy, hei_tangent = cos.get_splined_hei()

    # ts_geom = cos.images[0].copy()
    # ts_geom.coords = hei_coords
    # ts_geom.set_calculator(calc_getter())
    # print("Splined TS guess")
    # print(ts_geom.as_xyz())
    # with open("splined_ts_guess.xyz", "w") as handle:
        # handle.write(ts_geom.as_xyz())

    # print(f"Interpolated HEI: {hei_coords}")

    # ts_geom.set_calculator(calc_getter())
    # geoms = [ts_geom, ]
    # dimer_kwargs.update({
        # "restrict_step": "max",
        # "N_init": hei_tangent,
        # "calc_getter": calc_getter,
    # })
    # dimer_result = dimer_method(geoms, **dimer_kwargs)
    # dimer_cycles = dimer_result.dimer_cycles
    # last_cycle = dimer_cycles[-1]
    # ts_coords = last_cycle.trans_coords[1]
    # print(f"Optimized TS coord:")
    # print(ts_geom.as_xyz())


def run_calculations(geoms, calc_getter, path, calc_key, calc_kwargs,
                     scheduler=None, assert_track=False):
    print("Running calculations")
    def par_calc(geom):
        geom.calculator.run_calculation(geom.atoms, geom.coords)
        return geom

    for i, geom in enumerate(geoms):
        geom.set_calculator(calc_getter(i))
    if assert_track:
        assert all([geom.calculator.track for geom in geoms]), \
            "'track: True' must be present in calc section."

    if scheduler:
        client =  Client(scheduler, pure=False, silence_logs=False)
        geom_futures = client.map(par_calc, geoms)
        geoms = client.gather(geom_futures)
    else:
        for i, geom in enumerate(geoms):
            start = time.time()
            geom.calculator.run_calculation(geom.atoms, geom.coords)
            if i < (len(geoms)-2) and hasattr(geom.calculator, "propagate_wavefunction"):
                next_calculator = geoms[i+1].calculator
                geom.calculator.propagate_wavefunction(next_calculator)
            end = time.time()
            diff = end - start
            print(f"Ran calculation {i+1:02d}/{len(geoms):02d} in {diff:.1f} s.")
            sys.stdout.flush()
    return geoms


def get_overlapper(run_dict):
    try:
        calc_key = run_dict["calc"].pop("type")
    except KeyError:
        print("Creating Overlapper without calc_key.")
        calc_key = None
    calc_kwargs = run_dict["calc"]
    cwd = Path(".")
    ovlp_with = run_dict["overlaps"]["ovlp_with"]
    prev_n = run_dict["overlaps"]["prev_n"]
    overlapper = Overlapper(cwd, ovlp_with, prev_n, calc_key, calc_kwargs)
    return overlapper


def restore_calculators(run_dict):
    overlapper = get_overlapper(run_dict)

    cwd = Path(".")
    glob = run_dict["overlaps"]["glob"]
    regex = run_dict["overlaps"]["regex"]
    # First try globbing
    if glob:
        paths = natsorted([p for p in cwd.glob(glob)])
        if len(paths) == 0:
            raise Exception("Couldn't find any paths! Are you sure that your "
                           f"glob '{glob}' is right?")
        xyz_fns = [list(p.glob("*.xyz"))[0] for p in paths]
        geoms = [geom_from_xyz_file(xyz) for xyz in xyz_fns]
        if regex:
            mobjs = [re.search(regex, str(path)) for path in paths]
            calc_numbers = [int(mobj[1]) for mobj in mobjs]
            assert len(calc_numbers) == len(paths)
        else:
            calc_numbers = range(len(paths))
        [overlapper.set_files_from_dir(geom, path, calc_number)
         for calc_number, geom, path in zip(calc_numbers, geoms, paths)]
    else:
        # Otherwise check if geometries are defined in the run_dict
        if run_dict["xyz"]:
            geoms = get_geoms(run_dict["xyz"])
        else:
            # Else resort to globbing arbitrary xyz files
            xyz_fns = [str(p) for p in cwd.glob("*.xyz")]
            geoms = [geom_from_xyz_file(xyz) for xyz in xyz_fns]
        # geoms = geoms[:2]
        calc_num = overlapper.restore_calculators(geoms)
        geoms = geoms[:calc_num]
    return overlapper, geoms


def overlaps(run_dict, geoms=None):
    pickle_path = Path("pickles")
    if pickle_path.is_file() and confirm_input("Load pickled geoms?"):
        with open(pickle_path, "rb") as handle:
            overlapper, *geoms = cloudpickle.load(handle)
        print(f"Loaded overlap and {len(geoms)} from {str(pickle_path)}.")
    else:
        if not geoms:
            overlapper, geoms = restore_calculators(run_dict)
        else:
            overlapper = get_overlapper(run_dict)

        to_pickle = [overlapper] + geoms
        with open(pickle_path, "wb") as handle:
            cloudpickle.dump(to_pickle, handle)
    print("DEBUG Recreating overlapper!")
    overlapper = get_overlapper(run_dict)

    ovlp_dict = run_dict["overlaps"]
    ovlp_type = ovlp_dict["type"]
    double_mol = ovlp_dict["ao_ovlps"]
    recursive = ovlp_dict["recursive"]
    consider_first = ovlp_dict["consider_first"]
    skip = ovlp_dict["skip"]

    if ovlp_type == "wf" and double_mol:
        print("!"*10)
        print("WFOverlaps with true AO overlaps seem faulty right now!")
        print("!"*10)


    overlapper.overlaps_for_geoms(geoms,
                                  ovlp_type=ovlp_type,
                                  double_mol=double_mol,
                                  recursive=recursive,
                                  consider_first=consider_first,
                                  skip=skip,)

def run_opt(geom, calc_getter, opt_getter):
    geom.set_calculator(calc_getter(0))
    opt = opt_getter(geom)
    opt.run()
    return opt


def run_stocastic(stoc):#geom, calc_kwargs):
    # Fragment
    stoc.run()


def run_irc(geom, irc_kwargs):
    irc_type = irc_kwargs.pop("type")
    irc = IRC_DICT[irc_type](geom, **irc_kwargs)
    irc.run()


def run_tsopt(geom, tsopt_key, tsopt_kwargs):
    tsopt = TSOPT_DICT[tsopt_key](geom, **tsopt_kwargs)
    tsopt.run()


def copy_yaml_and_geometries(run_dict, yaml_fn, destination, new_yaml_fn=None):
    try:
        print(f"Trying to create directory '{destination}' ... ", end="")
        os.mkdir(destination)
        print("done")
    except FileExistsError:
        print("already exists")
    xyzs = run_dict["xyz"]
    print("Copying:")
    if isinstance(xyzs, str):
        xyzs = [xyzs, ]
    for xyz in xyzs:
        shutil.copy(xyz, destination)
        print("\t", xyz)
    shutil.copy(yaml_fn, destination)
    print("\t", yaml_fn)


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
            "type": None,
            "between": 0,
        },
        "cos": None,
        "dimer": None,
        "calc": {
            "pal": 1,
        },
        "opt": None,
        "tsopt": None,
        "overlaps": None,
        "glob": None,
        "stocastic": None,
        "xyz": None,
        "coord_type": "cart",
        "shake": None,
        "irc": None,
        "preopt_ends": False,
    }
    if "cos" in conf_dict:
        dd["cos"] = {
            "type": "neb",
            "parallel": 0,
            "fix_ends": True,
        }
        dd["opt"] = {
            "type": "qm",
            "align": True,
            "dump": True,
        }
    elif "opt" in conf_dict:
        dd["opt"] = {
            "type": "cg",
            "dump": True,
            "alpha": 0.25,
        }
    elif "overlaps" in conf_dict:
        dd["overlaps"] = {
            "type": "tden",
            "ao_ovlps": None,
            "glob": None,
            "recursive": False,
            "consider_first": None,
            "skip": 0,
            "regex": None,
            "ovlp_with": "previous",
            "prev_n": 0,
        }
    elif "stocastic" in conf_dict:
        dd["stocastic"] = {
            "type": "frag",
        }

    if "tsopt" in conf_dict:
        type_ = conf_dict["tsopt"]["type"]
        tsopt_dicts = {
            "dimer": {
                "type": "dimer",
                "max_step": 0.1,
                "max_cycles": 50,
                "trial_angle": 5,
                "angle_tol": 1,
                "dR_base": 0.01,
                "f_tran_mod": False,
                "multiple_translations": False,
            },
            "rsprfo": {
                "type": "rsprfo",
                # "max_step_length": .3,
                # "recalc_hess": None,
                "dump": True,
            },
        }
        tsopt_dict = tsopt_dicts[type_]
        tsopt_dict["do_hess"] = False
        dd["tsopt"] = tsopt_dict

    if "shake" in conf_dict:
        dd["shake"] = {
            "scale": 0.1,
            "seed": None,
        }

    if "irc" in conf_dict:
        dd["irc"] = {
            "type": "euler",
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
    yaml_dict = yaml.load(yaml_str, Loader=yaml.SafeLoader)
    # Load defaults to have a sane baseline
    run_dict = get_defaults(yaml_dict)
    # Update nested entries
    key_set = set(yaml_dict.keys())
    for key in key_set & set(("cos", "opt", "interpol", "overlaps",
                              "stocastic", "tsopt", "shake", "irc")):
        run_dict[key].update(yaml_dict[key])
    # Update non nested entries
    for key in key_set & set(("calc", "xyz", "pal", "coord_type",
                              "preopt_ends", )):
        run_dict[key] = yaml_dict[key]
    return run_dict


def dry_run(calc, geom):
    atoms, c3d = geom.atoms, geom.coords3d
    inp = calc.prepare_input(atoms, c3d.flatten(), "force")
    if not inp:
        print("Calculator does not use an explicit input file!")
        return
    with open(calc.inp_fn, "w") as handle:
        handle.write(inp)
    print(f"Wrote input to {calc.inp_fn}.")


def main(run_dict, restart=False, yaml_dir="./", scheduler=None,
         dryrun=None):
    xyz = run_dict["xyz"]
    if run_dict["interpol"]:
        interpolate = run_dict["interpol"]["type"]
        between = run_dict["interpol"]["between"]
    if run_dict["opt"]:
        opt_key = run_dict["opt"].pop("type")
        opt_kwargs = run_dict["opt"]
    if run_dict["cos"]:
        cos_key = run_dict["cos"].pop("type")
        cos_kwargs = run_dict["cos"]
        cos_kwargs["scheduler"] = scheduler
    if run_dict["stocastic"]:
        stoc_key = run_dict["stocastic"].pop("type")
        stoc_kwargs = run_dict["stocastic"]
    if run_dict["tsopt"]:
        tsopt_key = run_dict["tsopt"].pop("type")
        tsopt_kwargs = run_dict["tsopt"]
    if run_dict["irc"]:
        irc_kwargs = run_dict["irc"]

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

    # Prepare calculator
    calc_key = run_dict["calc"].pop("type")
    calc_kwargs = run_dict["calc"]
    calc_kwargs["out_dir"] = yaml_dir
    calc_getter = lambda index: get_calc(index, "image", calc_key, calc_kwargs)

    # Prepare optimizer
    opt_getter = lambda geoms: OPT_DICT[opt_key](geoms, **opt_kwargs)

    coord_type = run_dict["coord_type"]
    if run_dict["preopt_ends"]:
        # Update xyz list with optimized endpoint filenames
        xyz = preopt_ends(xyz, calc_getter)
    geoms = get_geoms(xyz, interpolate, between, coord_type=coord_type)
    if between and len(geoms) > 1:
        dump_geoms(geoms, "interpolated")

    if dryrun:
        calc = calc_getter(0)
        dry_run(calc, geoms[0])
        return
    elif run_dict["overlaps"]:
        geoms = run_calculations(geoms, calc_getter, yaml_dir, calc_key,
                                 calc_kwargs, scheduler, assert_track=True)
        overlaps(run_dict, geoms)
    elif run_dict["cos"]:
        cos_cls = COS_DICT[cos_key]
        if (issubclass(cos_cls, GrowingChainOfStates)
            or isinstance(cos_cls, type(FreezingString))):
            cos_kwargs["calc_getter"] = get_calc_closure("image", calc_key, calc_kwargs)
        cos = COS_DICT[cos_key](geoms, **cos_kwargs)
        run_cos(cos, calc_getter, opt_getter)
        if run_dict["tsopt"]:
            calc_getter = get_calc_closure(tsopt_key, calc_key, calc_kwargs)
            run_cos_tsopt(cos, tsopt_key, tsopt_kwargs, calc_getter)
    elif run_dict["opt"]:
        assert(len(geoms) == 1)
        geom = geoms[0]
        if run_dict["shake"]:
            shaked_coords = shake_coords(geom.coords, **run_dict["shake"])
            geom.coords = shaked_coords
        run_opt(geom, calc_getter, opt_getter)
    elif run_dict["stocastic"]:
        assert len(geoms) == 1
        geom = geoms[0]
        stoc_kwargs["calc_kwargs"] = calc_kwargs
        stoc = STOCASTIC_DICT[stoc_key](geom, **stoc_kwargs)
        run_stocastic(stoc)
    elif run_dict["irc"]:
        assert len(geoms) == 1
        geom = geoms[0]
        geom.set_calculator(calc_getter(0))
        run_irc(geom, irc_kwargs)
    elif run_dict["tsopt"]:
        assert len(geoms) == 1
        geom = geoms[0]
        geom.set_calculator(calc_getter(0))
        run_tsopt(geom, tsopt_key, tsopt_kwargs)
    else:
        geoms = run_calculations(geoms, calc_getter, yaml_dir, calc_key,
                                 calc_kwargs, scheduler)


def clean(force=False):
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
        "tsoptimizer.log",
        "wfoverlap.log",
        "host_*.calculator.log",
        "host_*.wfoverlap.log",
        "wfo_*.out"
        "optimization.trj",
        "cos.log",
        "*.gradient",
        "image_results.yaml",
        "optimizer_results.yaml",
        # ORCA specific
        "*orca.gbw",
        "*orca.cis",
        "*orca.engrad",
        "*orca.hessian",
        "*orca.inp",
        # OpenMOLCAS specific
        "image*.RasOrb",
        "image*.in",
        "image*.JobIph",
        "calculator*.out",
        "calculator*.JobIph",
        "calculator*.RasOrb",
        "*rasscf.molden",
        # Gaussian specific
        "image*.fchk",
        "image*.log",
        "image*.chk",
        "image*.com",
        "image*.*dump_635r",
        "image*ao_ovlp_rec"
        # Turbomole specific
        "image*.mos",
        "image*.alpha",
        "image*.beta",
        "image*.control",
        "image*.ciss_a",
        "calculator_*.control",
        "calculator_*.coord"
        "calculator_*.mos",
        "calculator_*.ciss_a",
        "image*.sing_a",
        "calculator*.sing_a",
        "*wavefunction.molden",
        "*input.xyz",
        "*.coord"
        # WFOverlap specific
        "wfo_*.*.out",
        # XTB specific
        "image*.grad",
        "calculator*.grad",
        "image_*.*.mos",
        "image*.*.ao_ovlp_rec",
        "splined_ts_guess.xyz",
        "dimer_ts.xyz",
        "dimer_pickle",
        "interpolated.geom_*.xyz",
        # Wavefunction overlap
        "wfo_*",
        "image*.molden",
        "crashed_*",
        "jmol.spt",
        # "overlap_data.h5",
        "*_CDD.png",
        "*_CDD.cub",
        "internal_coords.log",
        "hei_tangent",
        "optimization.trj",
        "splined_hei.xyz",
        "ts_opt.xyz",
        "final_geometry.xyz",
        "calculated_init_hessian",
    )
    to_rm_paths = list()
    for glob in rm_globs:
        to_rm_paths.extend(list(cwd.glob(glob)))
    to_rm_strs = [str(p) for p in to_rm_paths]
    for s in to_rm_strs:
        print(s)

    def delete():
        for p in to_rm_paths:
            try:
                os.remove(p)
            except OSError:
                # os.rmdir(p)
                shutil.rmtree(p)
            print(f"Deleted {p}")
    if force:
        delete()
        return
    # If we dont force the cleaning ask for confirmation first
    if confirm_input("Delete these files?"):
        delete()
    else:
        print("Aborting")


def print_header():
    """Generated from https://asciiartgen.now.sh/?s=pysisyphus&style=colossal"""
    logo = """                           d8b                            888
                           Y8P                            888
                                                          888
88888b.  888  888 .d8888b  888 .d8888b  888  888 88888b.  88888b.  888  888 .d8888b
888 "88b 888  888 88K      888 88K      888  888 888 "88b 888 "88b 888  888 88K
888  888 888  888 "Y8888b. 888 "Y8888b. 888  888 888  888 888  888 888  888 "Y8888b.
888 d88P Y88b 888      X88 888      X88 Y88b 888 888 d88P 888  888 Y88b 888      X88
88888P"   "Y88888  88888P' 888  88888P'  "Y88888 88888P"  888  888  "Y88888  88888P'
888           888                            888 888
888      Y8b d88P                       Y8b d88P 888
888       "Y88P"                         "Y88P"  888                            """
    version = f"Version {get_versions()['version']}"
    print(f"{logo}\n\n{version}\n")


def run():
    start_time = time.time()
    args = parse_args(sys.argv[1:])

    if args.yaml:
        yaml_dir = Path(os.path.abspath(args.yaml)).parent
        init_logging(yaml_dir, args.scheduler)
        with open(args.yaml) as handle:
            yaml_str = handle.read()
        run_dict = handle_yaml(yaml_str)
        sys.stdout.flush()

    if args.cp:
        copy_yaml_and_geometries(run_dict, args.yaml, args.cp)
        return
    elif args.clean:
        clean()
        return
    elif args.fclean:
        clean(force=True)
        return

    print_header()
    run_dict_without_none = {k: v for k, v in run_dict.items()
                             if v is not None}
    pprint(run_dict_without_none)
    print()

    if args.overlaps:
        overlaps(run_dict)
    elif args.couplings:
        couplings(args.couplings)
    elif args.sort_by_overlaps:
        energies_fn, max_ovlp_inds_fn = args.sort_by_overlaps
        consider_first = args.consider_first
        sort_by_overlaps(energies_fn, max_ovlp_inds_fn, consider_first)
    else:
        main(run_dict, args.restart, yaml_dir, args.scheduler, args.dryrun)
    end_time = time.time()
    duration = int(end_time - start_time)
    print(f"pysisyphus run took {duration}s.")

if __name__ == "__main__":
    run()
