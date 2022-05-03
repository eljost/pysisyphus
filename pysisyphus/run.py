import argparse
from collections import namedtuple
import copy
import datetime
import itertools as it
import os
from math import modf
from pathlib import Path
import platform
from pprint import pprint
import re
import shutil
import sys
import textwrap
import time

from distributed import Client
import numpy as np
import scipy as sp
import yaml

from pysisyphus import __version__
from pysisyphus.calculators import *
from pysisyphus.config import OUT_DIR_DEFAULT, p_DEFAULT, T_DEFAULT
from pysisyphus.cos import *
from pysisyphus.cos.GrowingChainOfStates import GrowingChainOfStates
from pysisyphus.color import bool_color
from pysisyphus.exceptions import HEIIsFirstOrLastException
from pysisyphus.dynamics import (
    get_mb_velocities_for_geom,
    mdp,
    md,
    get_colvar,
    Gaussian,
)
from pysisyphus.drivers import (
    relaxed_1d_scan,
    run_opt,
    run_precontr,
    run_perf,
    print_perf_results,
)
from pysisyphus.drivers.barriers import do_endopt_ts_barriers

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import (
    confirm_input,
    shake_coords,
    print_barrier,
    get_tangent_trj_str,
)
from pysisyphus.helpers_pure import (
    find_closest_sequence,
    merge_sets,
    recursive_update,
    highlight_text,
    approx_float,
    results_to_json,
)
from pysisyphus.intcoords import PrimitiveNotDefinedException
from pysisyphus.intcoords.setup import get_bond_mat
from pysisyphus.init_logging import init_logging
from pysisyphus.intcoords.PrimTypes import PrimTypes, normalize_prim_inputs
from pysisyphus.intcoords.helpers import form_coordinate_union
from pysisyphus.intcoords.setup import get_bond_sets
from pysisyphus.interpolate import interpolate_all
from pysisyphus.irc import *
from pysisyphus.io import save_hessian
from pysisyphus.stocastic import *
from pysisyphus.thermo import can_thermoanalysis
from pysisyphus.trj import get_geoms, dump_geoms, standardize_geoms
from pysisyphus.xyzloader import write_geoms_to_trj
from pysisyphus.yaml_mods import get_loader, UNITS


CALC_DICT = {
    "afir": AFIR,
    "composite": Composite,
    "conical": ConicalIntersection,
    "dftb+": DFTBp,
    "dimer": Dimer,
    "dummy": Dummy,
    "energymin": EnergyMin,
    "ext": ExternalPotential,
    "g09": Gaussian09,
    "g16": Gaussian16,
    "ipiserver": IPIServer,
    "mopac": MOPAC,
    "multi": MultiCalc,
    "oniom": ONIOM,
    "openmolcas": OpenMolcas,
    "orca": ORCA,
    "orca5": ORCA5,
    "psi4": Psi4,
    "pyxtb": PyXTB,
    "remote": Remote,
    "turbomole": Turbomole,
    "xtb": XTB,
}

try:
    from pysisyphus.calculators.PySCF import PySCF

    CALC_DICT["pyscf"] = PySCF
except ImportError:
    pass

try:
    from pysisyphus.calculators.QCEngine import QCEngine

    CALC_DICT["qcengine"] = QCEngine
except ImportError:
    pass

COS_DICT = {
    "neb": NEB.NEB,
    "aneb": AdaptiveNEB.AdaptiveNEB,
    "feneb": FreeEndNEB.FreeEndNEB,
    "szts": SimpleZTS.SimpleZTS,
    "gs": GrowingString.GrowingString,
    "fs": FreezingString.FreezingString,
}

IRC_DICT = {
    "dvv": DampedVelocityVerlet,
    "euler": Euler,
    "eulerpc": EulerPC,
    "gs": GonzalezSchlegel,
    "imk": IMKMod,
    "lqa": LQA,
    "modekill": ModeKill,
    "rk4": RK4,
}

STOCASTIC_DICT = {
    "frag": FragmentKick,
    "kick": Kick,
}


def parse_args(args):
    parser = argparse.ArgumentParser()

    action_group = parser.add_mutually_exclusive_group(required=True)
    action_group.add_argument(
        "yaml", nargs="?", help="Start pysisyphus with input from a YAML file."
    )
    action_group.add_argument(
        "--clean", action="store_true", help="Ask for confirmation before cleaning."
    )
    action_group.add_argument(
        "--fclean",
        action="store_true",
        help="Force cleaning without prior confirmation.",
    )
    action_group.add_argument(
        "--bibtex",
        action="store_true",
        help="Print bibtex string for pysisyphus paper.",
    )
    action_group.add_argument(
        "-v", "--version", action="store_true", help="Print pysisyphus version."
    )

    run_type_group = parser.add_mutually_exclusive_group(required=False)
    run_type_group.add_argument(
        "--restart",
        action="store_true",
        help="Continue a previously crashed/aborted/... pysisphus run.",
    )
    run_type_group.add_argument(
        "--cp",
        "--copy",
        nargs="+",
        help="Copy .yaml file and corresponding geometries from the 'geom' section "
        "to a new directory. The first argument is interpreted as destination. Any "
        "remaining (optional) arguments are files that are also copied.",
    )

    parser.add_argument(
        "--scheduler", default=None, help="Address of the dask scheduler."
    )
    return parser.parse_args()


def get_calc_closure(base_name, calc_key, calc_kwargs, iter_dict=None, index=None):
    if iter_dict is None:
        iter_dict = dict()

    if index is None:
        index = 0
    # Maps YAML input to actual Calculator-class arguments
    calc_map = {
        "calculator": "calculator",
        "calc": "calculator",  # shortcut for 'calculator'
        # calc1/calc2 are used for ConicalIntersection and EnergyMin.
        "calculator1": "calculator1",
        "calc1": "calculator1",  # shortcut for 'calculator1'
        "calcualtor2": "calculator2",
        "calc2": "calculator2",  # shortcut for 'calculator2'
    }

    def calc_getter(**add_kwargs):
        nonlocal index

        kwargs_copy = copy.deepcopy(calc_kwargs)

        # Some calculators are just wrappers, modifying forces from actual calculators,
        # e.g. AFIR and Dimer. If we find one of the keys in 'calc_map' in 'calc_kwargs'
        # we create the actual calculator and assign it to the corresponding value in
        # 'calc_map'.
        for key, val in calc_map.items():
            if key in kwargs_copy:
                # Use different base_name to distinguish the calculator(s)
                actual_base_name = val
                actual_kwargs = kwargs_copy.pop(key)
                actual_key = actual_kwargs.pop("type")
                # Pass 'index' to arguments, to avoid recreating calculators with
                # the same name.
                actual_calc = get_calc_closure(
                    actual_base_name, actual_key, actual_kwargs, index=index
                )()
                kwargs_copy[val] = actual_calc

        kwargs_copy["base_name"] = base_name
        kwargs_copy.update(add_kwargs)
        kwargs_copy["calc_number"] = index
        for key, iter_ in iter_dict.items():
            kwargs_copy[key] = next(iter_)
        index += 1
        return CALC_DICT[calc_key](**kwargs_copy)

    return calc_getter


def run_tsopt_from_cos(
    cos,
    tsopt_key,
    tsopt_kwargs,
    calc_getter=None,
    ovlp_thresh=0.4,
    coordinate_union="bonds",
):
    print(highlight_text(f"Running TS-optimization from COS"))

    # Later want a Cartesian HEI tangent, so if not already present we create
    # a Cartesian COS object to obtain the tangent from.
    atoms = cos.images[0].atoms
    if cos.coord_type != "cart":
        cart_images = list()
        for image in cos.images:
            cart_image = Geometry(atoms, image.cart_coords)
            cart_image.energy = image.energy
            cart_images.append(cart_image)
        cart_cos = ChainOfStates.ChainOfStates(cart_images)
    # Just continue using the Cartesian COS object
    else:
        cart_cos = cos

    hei_kind = tsopt_kwargs.pop("hei_kind", "splined")
    # Use plain, unsplined, HEI
    if hei_kind == "plain":
        hei_index = cos.get_hei_index()
        hei_image = cos.images[hei_index]
        # Select the Cartesian tangent from the COS
        cart_hei_tangent = cart_cos.get_tangent(hei_index)
    # Use splined HEI
    elif hei_kind == "splined":
        # The splined HEI tangent is usually very bady for the purpose of
        # selecting an imaginary mode to follow uphill. So we construct a better
        # HEI tangent by mixing the two tangents closest to the splined HEI.
        hei_coords, hei_energy, splined_hei_tangent, hei_index = cos.get_splined_hei()
        hei_image = Geometry(atoms, hei_coords)
        close_to_first = hei_index < 0.5
        close_to_last = hei_index > len(cos.images) - 1.5
        if close_to_first or close_to_last:
            close_to = "first" if close_to_first else "last"
            print(
                f"Splined HEI is too close to the {close_to} image. Aborting TS optimization!"
            )
            raise HEIIsFirstOrLastException
        # The hei_index is a float. We split off the decimal part and mix the two
        # nearest tangents accordingly.
        frac, floor = modf(hei_index)
        # Indices of the two nearest images with integer indices.
        floor = int(floor)
        ceil = floor + 1
        floor_tangent = cart_cos.get_tangent(floor)
        ceil_tangent = cart_cos.get_tangent(ceil)
        print(f"Creating mixed HEI tangent, using tangents at images {(floor, ceil)}.")
        print("Overlap of splined HEI tangent with these tangents:")
        for ind, tang in ((floor, floor_tangent), (ceil, ceil_tangent)):
            print(f"\t{ind:02d}: {splined_hei_tangent.dot(tang):.6f}")
        # When frac is big, e.g. 0.9 the tangent resembles the tangent at 'ceil'
        # so we mix in only (1-frac) == (1-0.9) == 0.1 of the 'floor' tangent.
        cart_hei_tangent = (1 - frac) * floor_tangent + frac * ceil_tangent
        cart_hei_tangent /= np.linalg.norm(cart_hei_tangent)
        # print(f"\t(1-{frac:.4f})*t({floor})+{frac:.4f}*t({ceil}): "
        # f"{splined_hei_tangent.dot(cart_hei_tangent):.6f}"
        # )
    else:
        raise Exception(f"Invalid hei_kind='{hei_kind}'!")

    print(f"Index of {hei_kind} highest energy image (HEI) is {hei_index:.2f}.")
    print()

    # When the COS was optimized in internal coordinates the united primitive
    # indices are already present and we just keep on using them.
    try:
        typed_prims = hei_image.internal.typed_prims
    # If the COS was optimized in Cartesians we have to generated a new
    # set of primitive internals.
    except AttributeError:

        def get_int_geom(geom):
            return Geometry(geom.atoms, geom.cart_coords, coord_type="redund")

        internal_geom1 = get_int_geom(cos.images[0])
        internal_geom2 = get_int_geom(cos.images[-1])
        typed_prims = form_coordinate_union(internal_geom1, internal_geom2)

    coord_kwargs = dict()
    if coordinate_union == "all":
        coord_kwargs["typed_prims"] = typed_prims
        union_msg = "Using full coordinate union for TS guess. Probably a bad idea!"
    elif coordinate_union in ("bonds", "bonds_bends_dihedrals"):
        # Only keep actual bonds ...
        valid_prim_types = (PrimTypes.BOND,)
        # ... and bends and dihedrals, if requested
        if coordinate_union == "bonds_bends_dihedrals":
            valid_prim_types += (PrimTypes.BEND, PrimTypes.PROPER_DIHEDRAL)
        coord_kwargs["define_prims"] = [
            tp for tp in typed_prims if tp[0] in valid_prim_types
        ]
        union_msg = f"Kept primitive types: {valid_prim_types}"
    else:
        union_msg = "No coordinate union."
    print(union_msg)

    ts_geom_kwargs = tsopt_kwargs.pop("geom")
    ts_coord_type = ts_geom_kwargs.pop("type")
    if ts_coord_type != "cart":
        ts_geom_kwargs["coord_kwargs"].update(coord_kwargs)

    ts_geom = Geometry(
        hei_image.atoms,
        hei_image.cart_coords,
        coord_type=ts_coord_type,
        **ts_geom_kwargs,
    )

    # Convert tangent from whatever coordinates to redundant internals.
    # When the HEI was splined the tangent will be in Cartesians.
    if ts_coord_type == "cart":
        ref_tangent = cart_hei_tangent
    elif ts_coord_type in ("redund", "dlc"):
        ref_tangent = ts_geom.internal.B @ cart_hei_tangent
    else:
        raise Exception(f"Invalid coord_type='{ts_coord_type}'!")
    ref_tangent /= np.linalg.norm(ref_tangent)

    # Dump HEI data
    #
    # Cartesian tangent and an animated .trj file
    cart_hei_fn = "cart_hei_tangent"
    np.savetxt(cart_hei_fn, cart_hei_tangent)
    trj = get_tangent_trj_str(
        ts_geom.atoms, ts_geom.cart_coords, cart_hei_tangent, points=10
    )
    trj_fn = cart_hei_fn + ".trj"
    with open(trj_fn, "w") as handle:
        handle.write(trj)
    print(f"Wrote animated HEI tangent to {trj_fn}\n")

    # Print HEI information (coords & tangent)
    print(f"{hei_kind.capitalize()} HEI (TS guess)")
    print(ts_geom.as_xyz())
    print()
    dummy = Geometry(atoms, cart_hei_tangent)
    print(f"{hei_kind.capitalize()} Cartesian HEI tangent")
    print(dummy.as_xyz())
    print()

    # Write out HEI information (coords & tangent)
    hei_xyz_fn = f"{hei_kind}_hei.xyz"
    with open(hei_xyz_fn, "w") as handle:
        handle.write(ts_geom.as_xyz())
    print(f"Wrote {hei_kind} HEI coordinates to '{hei_xyz_fn}'")

    ts_calc = calc_getter()

    def wrapped_calc_getter():
        return ts_calc

    ts_geom.set_calculator(ts_calc)
    tsopt_kwargs["prefix"] = "ts"

    if tsopt_key == "dimer":
        # Right now Dimer optimization is rectricted to cartesian
        # rotations and translations, even though translation in
        # internals would be possible.
        ts_geom = Geometry(hei_image.atoms, hei_image.cart_coords)
        dimer_kwargs = tsopt_kwargs.pop("dimer_kwargs", {})
        dimer_kwargs.update(
            {
                "N_raw": cart_hei_tangent,
                "base_name": "dimer",
            }
        )

        def wrapped_calc_getter():
            return Dimer(ts_calc, **dimer_kwargs)

        tsopt_key = "plbfgs"
        # When calling run_opt we pass the Hessian as additional argument,
        # so it is not recalculated unecessary. As no Hessian is available
        # for the dimer method we set it None.
        cart_hessian = None
    # Determine which imaginary mode has the highest overlap with the splined HEI tangent.
    else:
        print(f"Calculating Hessian at {hei_kind} TS guess.")
        # Calculate Hessian
        cart_hessian = ts_geom.cart_hessian
        # Continue Hessian in whatever coordinate system is actually in use
        H = ts_geom.hessian
        eigvals, eigvecs = np.linalg.eigh(H)
        neg_inds = eigvals < -1e-4
        if sum(neg_inds) == 0:
            raise Exception("No negative eigenvalues found at splined HEI. Exiting!")
        eigval_str = np.array2string(eigvals[neg_inds], precision=6)
        print(f"Negative eigenvalues at splined HEI:\n{eigval_str}")
        neg_eigvals = eigvals[neg_inds]
        neg_eigvecs = eigvecs.T[neg_inds]
        ovlps = [np.abs(imag_mode.dot(ref_tangent)) for imag_mode in neg_eigvecs]
        print("Overlaps between HEI tangent and imaginary modes:")
        for i, ov in enumerate(ovlps):
            print(f"\t{i:02d}: {ov:.6f}")
        max_ovlp_ind = np.argmax(ovlps)
        max_ovlp = ovlps[max_ovlp_ind]
        print(
            f"Imaginary mode {max_ovlp_ind} has highest overlap ({max_ovlp:.2%}) "
            "with splined HEI tangent."
        )
        rel_ovlps = np.array(ovlps) / max(ovlps)
        similar_inds = rel_ovlps > 0.80
        # Only 1 big overlap is present
        if (max_ovlp >= ovlp_thresh) and (similar_inds.sum() == 1):
            ovlp_root = np.argmax(ovlps)
        # Multiple big and similar overlaps are present.
        elif (max_ovlp >= ovlp_thresh) and (similar_inds.sum() > 1):
            # Will yield the first occurence of True, which corresponds to a
            # similar overlaps with the most negative eigenvalue.
            ovlp_root = similar_inds.argmax()
            neg_eigval = neg_eigvals[ovlp_root]
            verbose_inds = np.arange(neg_eigvals.size)[similar_inds]
            print(
                f"Overlaps {verbose_inds} are very similar! Falling back to the "
                f"one with the most negative eigenvalue {neg_eigval:.6f} "
                f"(mode {ovlp_root})."
            )
        # Fallback to the most negative eigenvalue when all overlaps are too low.
        else:
            ovlp_root = neg_eigvals.argmin()
            neg_eigval = neg_eigvals[ovlp_root]
            print(
                f"Highest overlap {max_ovlp:.6f} is below the threshold "
                f"of {ovlp_thresh:.6f}.\nFalling back to mode {ovlp_root} with most "
                f"negative eigenvalue {neg_eigval:.6f}."
            )
        root = tsopt_kwargs.get("root", None)
        if root is None:
            # Use mode with highest overlap as initial root
            tsopt_kwargs["root"] = ovlp_root
        else:
            print(
                f"Initial root={root} given, neglecting root {ovlp_root} "
                "determined from overlaps."
            )

    opt_result = run_opt(
        ts_geom,
        calc_getter=wrapped_calc_getter,
        opt_key=tsopt_key,
        opt_kwargs=tsopt_kwargs,
        cart_hessian=cart_hessian,
        title="TS-Optimization",
        copy_final_geom="ts_opt.xyz",
    )
    ts_geom = opt_result.geom
    ts_opt = opt_result.opt

    # Restore original calculator for Dimer calculations
    if tsopt_key == "dimer":
        ts_geom.set_calculator(ts_calc)

    print(f"Optimized TS coords:")
    print(ts_geom.as_xyz())
    # Include ts_ prefix
    ts_opt_fn = ts_opt.get_path_for_fn("opt.xyz")
    with open(ts_opt_fn, "w") as handle:
        handle.write(ts_geom.as_xyz())
    print(f"Wrote TS geometry to '{ts_opt_fn}\n'")

    ts_energy = ts_geom.energy
    first_cos_energy = cos.images[0].energy
    last_cos_energy = cos.images[-1].energy
    print_barrier(ts_energy, first_cos_energy, "TS", "first COS image")
    print_barrier(ts_energy, last_cos_energy, "TS", "last COS image")

    print()
    return opt_result


def run_calculations(
    geoms,
    calc_getter,
    scheduler=None,
    assert_track=False,
    run_func=None,
):
    print(highlight_text("Running calculations"))

    func_name = "run_calculation" if run_func is None else run_func

    def par_calc(geom):
        return getattr(geom.calculator, func_name)(geom.atoms, geom.coords)

    for geom in geoms:
        geom.set_calculator(calc_getter())

    if assert_track:
        assert all(
            [geom.calculator.track for geom in geoms]
        ), "'track: True' must be present in calc section."

    if scheduler:
        client = Client(scheduler, pure=False, silence_logs=False)
        results_futures = client.map(par_calc, geoms)
        all_results = client.gather(results_futures)
    else:
        all_results = list()
        i_fmt = "02d"
        for i, geom in enumerate(geoms):
            print(highlight_text(f"Calculation {i:{i_fmt}}", level=1))

            start = time.time()
            print(geom)
            results = getattr(geom.calculator, func_name)(geom.atoms, geom.cart_coords)
            # results dict of MultiCalc will contain keys that can be dumped yet. So
            # we skip the JSON dumping when KeyError is raised.
            try:
                as_json = results_to_json(results)
                calc = geom.calculator
                # Decrease counter, because it will be increased by 1, w.r.t to the
                # calculation.
                json_fn = calc.make_fn("results", counter=calc.calc_counter - 1)
                with open(json_fn, "w") as handle:
                    handle.write(as_json)
            except KeyError:
                print("Skipped JSON dump of calculation results!")

            hess_keys = [
                key
                for key, val in results.items()
                if isinstance(val, dict) and "hessian" in val
            ]
            for hkey in hess_keys:
                hres = results[hkey]
                hfn = f"{hkey}_hessian.h5"
                save_hessian(
                    hfn,
                    geom,
                    cart_hessian=hres["hessian"],
                    energy=hres["energy"],
                )
                print(f"Dumped hessian to '{hfn}'.")

            all_results.append(results)
            if i < (len(geoms) - 1):
                try:
                    cur_calculator = geom.calculator
                    next_calculator = geoms[i + 1].calculator
                    next_calculator.set_chkfiles(cur_calculator.get_chkfiles())
                    msg = f"Set chkfiles on calculator {i:{i_fmt}}"
                except AttributeError:
                    msg = "Calculator does not support set/get_chkfiles!"
                print(msg)
            end = time.time()
            diff = end - start
            print(f"Calculation took {diff:.1f} s.\n")
            sys.stdout.flush()
    print()

    for geom, results in zip(geoms, all_results):
        try:
            geom.set_results(results)
        except KeyError:
            pass

    return geoms, all_results


def run_stocastic(stoc):
    # Fragment
    stoc.run()
    print()

    return stoc


def run_md(geom, calc_getter, md_kwargs):
    print(highlight_text(f"Running Molecular Dynamics"))

    calc = calc_getter()
    geom.set_calculator(calc)

    T_init_vel = md_kwargs.pop("T_init_vel")
    steps = md_kwargs.pop("steps")
    dt = md_kwargs.pop("dt")
    seed = md_kwargs.pop("seed", None)

    _gaussian = md_kwargs.pop("gaussian", {})
    gaussians = list()
    for g_name, g_kwargs in _gaussian.items():
        # Create collective variable
        cv_kwargs = g_kwargs.pop("colvar")
        cv_key = cv_kwargs.pop("type")
        colvar = get_colvar(cv_key, cv_kwargs)

        # Create & append Gaussian
        g_w = g_kwargs.pop("w")
        g_s = g_kwargs.pop("s")
        g_stride = g_kwargs.pop("stride")
        gau = Gaussian(w=g_w, s=g_s, colvar=colvar, dump_name=g_name)
        gaussians.append((g_name, gau, g_stride))

    remove_com_v = md_kwargs.get("remove_com_v")
    v0 = get_mb_velocities_for_geom(
        geom, T_init_vel, seed=seed, remove_com_v=remove_com_v, remove_rot_v=False
    ).flatten()
    md_result = md(geom, v0=v0, steps=steps, dt=dt, gaussians=gaussians, **md_kwargs)

    # from pysisyphus.xyzloader import coords_to_trj
    # trj_fn = "md.trj"
    # _ = coords_to_trj(
    # trj_fn, geom.atoms, md_result.coords[::md_kwargs["dump_stride"]]
    # )
    print()

    return md_result


def run_scan(geom, calc_getter, scan_kwargs, callback=None):
    print(highlight_text("Relaxed Scan") + "\n")
    assert (
        geom.coord_type != "cart"
    ), "Internal coordinates are required for coordinate scans."

    type_ = scan_kwargs["type"]
    indices = scan_kwargs["indices"]
    steps = scan_kwargs["steps"]
    start = scan_kwargs.get("start", None)
    end = scan_kwargs.get("end", None)
    step_size = scan_kwargs.get("step_size", None)
    symmetric = scan_kwargs.get("symmetric", False)
    # The final prim value is determined either as
    #  start + steps*step_size
    # or
    #  (end - start) / steps .
    #
    # So we always require steps and either end or step_size.
    # bool(a) != bool(b) amounts to an logical XOR.
    assert (steps > 0) and (
        bool(end) != bool(step_size)
    ), "Please specify either 'end' or 'step_size'!"
    if symmetric:
        assert step_size and (
            start is None
        ), "'symmetric: True' requires 'step_size' and 'start == None'!"

    constrain_prims = normalize_prim_inputs(((type_, *indices),))
    constr_prim = constrain_prims[0]

    start_was_none = start is None
    if start_was_none:
        try:
            constr_ind = geom.internal.get_index_of_typed_prim(constr_prim)
        # Recreate geom with appropriate primitive
        except PrimitiveNotDefinedException:
            geom = geom.copy(coord_kwargs={"define_prims": constrain_prims})
            constr_ind = geom.internal.get_index_of_typed_prim(constr_prim)
        # The given indices may not correspond exactly to a typed primitives,
        # as they may be reversed. So we fetch the actual typed primitive.
        constr_prim = geom.internal.typed_prims[constr_ind]
        start = geom.coords[constr_ind]

    if step_size is None:
        step_size = (end - start) / steps
    opt_kwargs = scan_kwargs["opt"].copy()
    opt_key = opt_kwargs.pop("type")

    def wrapper(geom, start, step_size, steps, pref=None):
        return relaxed_1d_scan(
            geom,
            calc_getter,
            [
                constr_prim,
            ],
            start,
            step_size,
            steps,
            opt_key,
            opt_kwargs,
            pref=pref,
            callback=callback,
        )

    if not symmetric:
        scan_geoms, scan_vals, scan_energies = wrapper(geom, start, step_size, steps)
    else:
        # Negative direction
        print(highlight_text("Negative direction", level=1) + "\n")
        minus_geoms, minus_vals, minus_energies = wrapper(
            geom, start, -step_size, steps, pref="minus"
        )
        init_geom = minus_geoms[0].copy()
        # Positive direction. Compared to the negative direction we start at a
        # displaced geometry and reduce the number of steps by 1.
        print(highlight_text("Positive direction", level=1) + "\n")
        plus_start = start + step_size
        # Do one step less, as we already start from the optimized geometry
        plus_steps = steps - 1
        plus_geoms, plus_vals, plus_energies = wrapper(
            init_geom, plus_start, step_size, plus_steps, pref="plus"
        )
        scan_geoms = minus_geoms[::-1] + plus_geoms
        scan_vals = np.concatenate((minus_vals[::-1], plus_vals))
        scan_energies = np.concatenate((minus_energies[::-1], plus_energies))

        trj = "\n".join([geom.as_xyz() for geom in scan_geoms])
        with open("relaxed_scan.trj", "w") as handle:
            handle.write(trj)
        scan_data = np.stack((scan_vals, scan_energies), axis=1)
        np.savetxt(f"relaxed_scan.dat", scan_data)
    return scan_geoms, scan_vals, scan_energies


def run_preopt(first_geom, last_geom, calc_getter, preopt_key, preopt_kwargs):
    strict = preopt_kwargs.pop("strict", False)

    geom_kwargs = preopt_kwargs.pop("geom")
    coord_type = geom_kwargs.pop("type")

    def recreate_geom(geom):
        return Geometry(geom.atoms, geom.coords, coord_type=coord_type, **geom_kwargs)

    first_geom = recreate_geom(first_geom)
    last_geom = recreate_geom(last_geom)

    opt_results = list()
    for geom, key in ((first_geom, "first"), (last_geom, "last")):
        # Backup original geometry for RMSD calculation
        org_geom = geom.copy(coord_type="cart")

        prefix = f"{key}_pre"
        opt_kwargs = preopt_kwargs.copy()
        opt_kwargs.update(
            {
                "prefix": prefix,
                "h5_group_name": prefix,
            }
        )
        opt_result = run_opt(
            geom, calc_getter, preopt_key, opt_kwargs, title=f"{key} preoptimization"
        )
        opt = opt_result.opt
        opt_results.append(opt_result)

        # Continue with next pre-optimization when stopped manually
        if strict and not opt.stopped and not opt.is_converged:
            print(f"Problem in preoptimization of {key}. Exiting!")
            sys.exit()
        print(f"Preoptimization of {key} geometry converged!")
        opt_fn = opt.get_path_for_fn(f"opt.xyz")
        shutil.move(opt.final_fn, opt_fn)
        print(f"Saved final preoptimized structure to '{opt_fn}'.")

        rmsd = org_geom.rmsd(opt_result.geom)
        print(f"RMSD with initial geometry: {rmsd:.6f} au")
        print()

    return opt_results


def run_irc(geom, irc_key, irc_kwargs, calc_getter):
    print(highlight_text(f"Running IRC") + "\n")

    calc = calc_getter()
    calc.base_name = "irc"
    geom.set_calculator(calc, clear=False)

    # Recreate geometry with Cartesian coordinates if needed.
    if geom.coord_type != "cart":
        geom = geom.copy_all(coord_type="cart")

    irc = IRC_DICT[irc_key](geom, **irc_kwargs)
    irc.run()

    return irc


def do_rmsds(xyz, geoms, end_fns, end_geoms, preopt_map=None, similar_thresh=0.25):
    if len(end_fns) != 2 or len(end_geoms) != 2:
        return
    end_fns = [str(fn) for fn in end_fns]
    max_end_len = max(len(s) for s in end_fns)

    if len(geoms) == 1:
        return
    elif len(geoms) > 2:
        geoms = (geoms[0], geoms[-1])
    assert len(geoms) == 2

    if isinstance(xyz, str):
        xyz = (f"{xyz}, first entry", f"{xyz}, last entry")
    elif (not isinstance(xyz, str)) and len(xyz) >= 2:
        xyz = (xyz[0], xyz[-1])
    assert len(xyz) == 2

    # Append original filenames, if supplied
    if preopt_map is not None:
        xyz = [f"{preopt}/{preopt_map[preopt]}" for preopt in xyz]

    max_len = max(len(s) for s in xyz)

    print(highlight_text(f"RMSDs After End Optimizations"))
    print()

    start_cbms = [get_bond_mat(geom) for geom in geoms]
    end_cbms = [get_bond_mat(geom) for geom in end_geoms]
    for i, start_geom in enumerate(geoms):
        fn = xyz[i]
        found_similar = False
        print(f"start geom {i:>2d} ({fn:>{max_len}s})")
        # Condensed bond mat
        start_cbm = start_cbms[i]
        for j, end_geom in enumerate(end_geoms):
            end_fn = end_fns[j]
            # Compare bond matrices
            cbm_match = (start_cbm == end_cbms[j]).all()
            cbm_str = "(bond matrices match)" if cbm_match else ""
            rmsd = start_geom.rmsd(end_geom)
            similar_str = ""
            if rmsd < similar_thresh:
                found_similar = True
                similar_str = " (similar)"
            print(
                f"\tend geom {j:>2d} ({end_fn:>{max_end_len}s}): "
                f"RMSD={rmsd:>8.6f} au{similar_str} " + cbm_str
            )
        if not found_similar:
            print(f"\tRMSDs of end geometries are dissimilar to '{fn}'!")
    print()


def run_endopt(irc, endopt_key, endopt_kwargs, calc_getter):
    print(highlight_text(f"Optimizing reaction path ends"))

    # Gather geometries that shall be optimized and appropriate keys.
    to_opt = list()
    if irc.forward:
        coords = irc.all_coords[0]
        to_opt.append((coords, "forward"))
    if irc.backward:
        coords = irc.all_coords[-1]
        to_opt.append((coords, "backward"))
    if irc.downhill:
        coords = irc.all_coords[-1]
        to_opt.append((coords, "downhill"))

    def to_frozensets(sets):
        return [frozenset(_) for _ in sets]

    separate_fragments = endopt_kwargs.pop("fragments", False)
    total = separate_fragments in ("total", False)

    # Convert to array for easy indexing with the fragment lists
    atoms = np.array(irc.atoms)
    fragments_to_opt = list()
    # Expand endpoints into fragments if requested
    for coords, key in to_opt:
        base_name = f"{key}_end"
        c3d = coords.reshape(-1, 3)
        fragments = list()
        fragment_names = list()

        # Detect separate fragments if requested.
        if separate_fragments:
            bond_sets = to_frozensets(get_bond_sets(atoms.tolist(), c3d))
            # Sort atom indices, so the atoms don't become totally scrambled.
            fragments.extend([sorted(frag) for frag in merge_sets(bond_sets)])
            # Disable higher fragment counts. I'm looking forward to the day
            # this ever occurs and someone complains :)
            assert len(fragments) < 10, "Something probably went wrong"
            fragment_names.extend(
                [f"{base_name}_frag{i:02d}" for i, _ in enumerate(fragments)]
            )
            print(f"Found {len(fragments)} fragment(s) at {base_name}")
            for frag_name, frag in zip(fragment_names, fragments):
                print(f"\t{frag_name}: {len(frag)} atoms")

        # Optimize the whole geometries, without splitting them into fragments
        # Skip this optimization if separate fragments are requested and only
        # one fragment is present, which would result in twice the same optimization.
        skip_one_frag = separate_fragments and len(fragments) == 1
        if total and not skip_one_frag:
            # Dummy fragment containing all atom indices.
            fragments.extend(
                [
                    range(len(atoms)),
                ]
            )
            fragment_names.extend(
                [
                    base_name,
                ]
            )
        elif skip_one_frag:
            print(
                f"Only one fragment present for '{key}'. Skipping optimization "
                "of total system."
            )

        fragment_keys = [key] * len(fragments)
        fragment_atoms = [tuple(atoms[list(frag)]) for frag in fragments]
        fragment_coords = [c3d[frag].flatten() for frag in fragments]
        fragments_to_opt.extend(
            list(zip(fragment_keys, fragment_names, fragment_atoms, fragment_coords))
        )
        print()

    to_opt = fragments_to_opt

    geom_kwargs = endopt_kwargs.pop("geom")
    coord_type = geom_kwargs.pop("type")

    results = {k: list() for k in ("forward", "backward", "downhill")}
    for key, name, atoms, coords in to_opt:
        geom = Geometry(
            atoms,
            coords,
            coord_type=coord_type,
            **geom_kwargs,
        )
        initial_fn = f"{name}_initial.xyz"
        with open(initial_fn, "w") as handle:
            handle.write(geom.as_xyz())

        def wrapped_calc_getter():
            calc = calc_getter()
            calc.base_name = name
            return calc

        opt_kwargs = endopt_kwargs.copy()
        opt_kwargs.update(
            {
                "prefix": name,
                "h5_group_name": name,
                "dump": True,
            }
        )
        try:
            opt_result = run_opt(
                geom,
                wrapped_calc_getter,
                endopt_key,
                opt_kwargs,
                title=f"{name} Optimization",
                level=1,
            )
        except Exception as err:
            print(f"{err}\nOptimization crashed!")
            continue
        final_fn = opt_result.opt.final_fn
        opt_fn = f"{name}_opt.xyz"
        shutil.move(final_fn, opt_fn)
        print(f"Moved '{final_fn.name}' to '{opt_fn}'.\n")
        results[key].append(opt_result)
    print()
    return results["forward"], results["backward"], results["downhill"]


def run_mdp(geom, calc_getter, mdp_kwargs):
    cwd = Path(".").resolve()
    for i in range(3):
        out_dir = cwd / f"outdir{i:02d}"
        try:
            os.mkdir(out_dir)
        except FileExistsError:
            print(out_dir, "exists")
        add_kwargs = {
            "out_dir": out_dir,
        }
        calc = calc_getter(**add_kwargs)
        geom.set_calculator(calc)
        mdp_result = mdp(geom, **mdp_kwargs)
        break
    return mdp_result


def copy_yaml_and_geometries(run_dict, yaml_fn, dest_and_add_cp, new_yaml_fn=None):
    destination, *copy_also = dest_and_add_cp
    src_path = Path(yaml_fn).resolve().parent
    destination = Path(destination)
    try:
        print(f"Trying to create directory '{destination}' ... ", end="")
        os.mkdir(destination)
        print("done")
    except FileExistsError:
        print("already exists")
    if "geom" in run_dict:
        xyzs = run_dict["geom"]["fn"]
    else:
        xyzs = run_dict["xyz"]
    print("Copying:")
    # Copy geometries
    # When newlines are present we have an inline xyz formatted string
    if not "\n" in xyzs:
        if isinstance(xyzs, str):
            xyzs = [
                xyzs,
            ]
        for xyz in xyzs:
            if xyz.startswith("lib:"):
                continue
            shutil.copy(src_path / xyz, destination)
            print("\t", xyz)
    else:
        print("Found inline xyz formatted string. No files to copy!")
    # Copy additional files
    for src in copy_also:
        try:
            shutil.copy(src_path / src, destination)
            print(f"\t{src}")
        except FileNotFoundError:
            print(f"\tCould not find '{src}'. Skipping!")
    # Update yaml_fn to match destination
    yaml_dest_fn = Path(destination.stem).with_suffix(".yaml")
    shutil.copy(yaml_fn, destination / yaml_dest_fn)
    print("\t", yaml_fn)


def get_defaults(conf_dict, T_default=T_DEFAULT, p_default=p_DEFAULT):
    # Defaults
    dd = {
        "assert": None,
        "barrier": None,
        "calc": {
            "type": "dummy",
        },
        "cos": None,
        "endopt": None,
        "geom": None,
        "glob": None,
        "interpol": None,
        "irc": None,
        "md": None,
        "mdp": None,
        "opt": None,
        "perf": None,
        "precontr": None,
        "preopt": None,
        "scan": None,
        "stocastic": None,
        "shake": None,
        "tsopt": None,
    }

    mol_opt_defaults = {
        "dump": True,
        "max_cycles": 150,
        "overachieve_factor": 5,
        "type": "rfo",
        "do_hess": False,
        "T": T_default,
        "p": p_default,
    }
    cos_opt_defaults = {
        "type": "qm",
        "align": True,
        "dump": True,
    }
    if "interpol" in conf_dict:
        dd["interpol"] = {
            "align": True,
        }

    if "cos" in conf_dict:
        dd["cos"] = {
            "type": "neb",
            "fix_ends": True,
        }
        dd["opt"] = cos_opt_defaults.copy()
    # Use a different, more powerful, optimizer when we are not dealing
    # with a COS-optimization.
    elif "opt" in conf_dict:
        dd["opt"] = mol_opt_defaults.copy()
    elif "stocastic" in conf_dict:
        dd["stocastic"] = {
            "type": "frag",
        }

    def get_opt_geom_defaults():
        return {
            "type": "redund",
            "coord_kwargs": {},
        }

    if "tsopt" in conf_dict:
        dd["tsopt"] = mol_opt_defaults.copy()
        dd["tsopt"].update(
            {
                "type": "rsprfo",
                "h5_group_name": "tsopt",
                "prefix": "ts",
            }
        )
        if "cos" in conf_dict:
            dd["tsopt"]["geom"] = get_opt_geom_defaults()

    if "precontr" in conf_dict:
        dd["precontr"] = {
            "prefix": "precontr",
        }

    if "preopt" in conf_dict:
        # We can't just copy dd["opt"] because there will probably be
        # some COS specific optimizer, but we just wan't to optimize the
        # (molecular) endpoints.
        dd["preopt"] = mol_opt_defaults.copy()
        dd["preopt"].update(
            {
                # Optimization specific
                # We are a bit more cautious here
                "max_cycles": 100,
                "thresh": "gau_loose",
                "trust_max": 0.3,
                "geom": get_opt_geom_defaults(),
                # Preopt specific
                "strict": False,
            }
        )
        # dd["preopt"]["geom"]["type"] = "tric"

    if "endopt" in conf_dict:
        dd["endopt"] = mol_opt_defaults.copy()
        dd["endopt"].update(
            {
                "thresh": "gau",
                "fragments": False,
                "geom": get_opt_geom_defaults(),
            }
        )

    if "barriers" in conf_dict:
        dd["barriers"] = {
            "T": T_default,
            "p": p_default,
            "solv_calc": {},
            "do_standard_state_corr": True,
        }

    if "shake" in conf_dict:
        dd["shake"] = {
            "scale": 0.1,
            "seed": None,
        }

    if "irc" in conf_dict:
        dd["irc"] = {
            "type": "eulerpc",
            "rms_grad_thresh": 1e-3,
        }

    if "assert" in conf_dict:
        dd["assert"] = {}

    if "geom" in conf_dict:
        dd["geom"] = {
            "type": "cart",
            "coord_kwargs": {},
        }

    if "mdp" in conf_dict:
        dd["mdp"] = {}

    if "scan" in conf_dict:
        dd["scan"] = {
            "opt": mol_opt_defaults.copy(),
            "symmetric": False,
        }
        dd["scan"]["opt"]["dump"] = False

    if "md" in conf_dict:
        md_T = T_default
        dd["md"] = {
            "T": md_T,
            "T_init_vel": md_T,
            "dt": 0.5,
            "thermostat": "csvr_2",
            "timecon": 50,
            "print_stride": 100,
            "dump_stride": 10,
            "remove_com_v": True,
        }

    if "perf" in conf_dict:
        dd["perf"] = {
            "mems": 2500,
            "repeat": 3,
        }

    return dd


def get_last_calc_cycle():
    def keyfunc(path):
        return re.match(r"image_\d+.(\d+).out", str(path))[1]

    cwd = Path(".")
    calc_logs = [str(cl) for cl in cwd.glob("image_*.*.out")]
    calc_logs = sorted(calc_logs, key=keyfunc)
    grouped = it.groupby(calc_logs, key=keyfunc)
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


def setup_run_dict(run_dict):
    org_dict = run_dict.copy()

    # Load defaults to have a sane baseline
    run_dict = get_defaults(run_dict)
    # Update nested entries that are dicts by themselves
    # Take care to insert a , after the string!
    key_set = set(org_dict.keys())
    for key in key_set & set(
        (
            "assert",
            "barriers",
            "calc",
            "cos",
            "endopt",
            "geom",
            "interpol",
            "irc",
            "md",
            "mdp",
            "opt",
            "perf",
            "precontr",
            "preopt",
            "scan",
            "shake",
            "stocastic",
            "tsopt",
        )
    ):
        try:
            # Recursive update, because there may be nested dicts
            recursive_update(run_dict[key], org_dict[key])
        except TypeError:
            print(f"Using default values for '{key}' section.")
    return run_dict


RunResult = namedtuple(
    "RunResult",
    (
        "preopt_first_geom preopt_last_geom "
        "cos cos_opt "
        "ts_geom ts_opt "
        "end_geoms irc irc_geom "
        "mdp_result "
        "opt_geom opt "
        "calced_geoms calced_results "
        "stocastic calc_getter "
        "scan_geoms scan_vals scan_energies "
        "perf_results "
    ),
)


def main(run_dict, restart=False, yaml_dir="./", scheduler=None):

    # Dump run_dict
    run_dict_copy = run_dict.copy()
    run_dict_copy["version"] = __version__
    with open("RUN.yaml", "w") as handle:
        yaml.dump(run_dict_copy, handle)

    if run_dict["interpol"]:
        interpol_key = run_dict["interpol"].pop("type")
        interpol_kwargs = run_dict["interpol"]
    # Preoptimization prior to COS optimization
    if run_dict["preopt"]:
        preopt_key = run_dict["preopt"].pop("type")
        preopt_kwargs = run_dict["preopt"]
    # Optimization of fragments after IRC integration
    if run_dict["endopt"]:
        endopt_key = run_dict["endopt"].pop("type")
        endopt_kwargs = run_dict["endopt"]
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
        irc_key = run_dict["irc"].pop("type")
        irc_kwargs = run_dict["irc"]

    # Handle geometry input. This section must always be present.
    geom_kwargs = run_dict["geom"]
    xyz = geom_kwargs.pop("fn")
    coord_type = geom_kwargs.pop("type")
    union = geom_kwargs.pop("union", None)

    ####################
    # CALCULATOR SETUP #
    ####################

    # Prepare calculator
    calc_key = run_dict["calc"].pop("type")
    calc_kwargs = run_dict["calc"]
    calc_run_func = calc_kwargs.pop("run_func", None)
    calc_kwargs["out_dir"] = calc_kwargs.get("out_dir", yaml_dir / OUT_DIR_DEFAULT)
    if calc_key in ("oniom", "ext"):
        geoms = get_geoms(xyz, quiet=True)
        iter_dict = {
            "geom": iter(geoms),
        }
    elif calc_key == "multi":
        geoms = get_geoms(xyz, quiet=True)
        iter_dict = {
            "base_name": iter([geom.name for geom in geoms]),
        }
    else:
        iter_dict = None
    calc_getter = get_calc_closure(
        "calculator", calc_key, calc_kwargs, iter_dict=iter_dict
    )
    # Create second function that returns a wrapped calculator. This may be
    # useful if we later want to drop the wrapper and use the actual calculator.
    if "calc" in calc_kwargs:
        act_calc_kwargs = calc_kwargs["calc"].copy()
        act_calc_key = act_calc_kwargs.pop("type")
        act_calc_getter = get_calc_closure(
            "act_calculator", act_calc_key, act_calc_kwargs
        )
    try:
        solv_calc_kwargs = run_dict["barriers"].pop("solv_calc")
        solv_calc_key = solv_calc_kwargs.pop("type")
        solv_calc_getter = get_calc_closure(
            "solv_calculator", solv_calc_key, solv_calc_kwargs
        )
    except KeyError:
        solv_calc_getter = None

    ##################
    # GEOMETRY SETUP #
    ##################

    # Initial loading of geometries from file(s)
    geoms = get_geoms(xyz, coord_type="cart")

    # ------------------------+
    #   Preconditioning of   |
    # Translation & Rotation |
    # ------------------------+

    if run_dict["precontr"]:
        ptr_geom0, ptr_geom_m1 = run_precontr(
            geoms[0], geoms[1], prefix=run_dict["precontr"]["prefix"]
        )
        geoms[0] = ptr_geom0
        geoms[-1] = ptr_geom_m1

    # -----------------------+
    # Preoptimization of    |
    # first & last image(s) |
    # -----------------------+

    # Preoptimization only makes sense with a subsequent COS run.
    if run_dict["preopt"] and run_dict["cos"]:
        assert len(geoms) > 1
        first_opt_result, last_opt_result = run_preopt(
            geoms[0],
            geoms[-1],
            calc_getter,
            preopt_key,
            preopt_kwargs,
        )
        # Update with (pre)optimized geometries. 'preopt_first_geom'/'preopt_last_geom'
        # are assigned here so they can later be assigned to 'run_result':
        geoms[0] = preopt_first_geom = first_opt_result.geom.copy(coord_type=coord_type)
        geoms[-1] = preopt_last_geom = last_opt_result.geom.copy(coord_type=coord_type)

    if run_dict["interpol"]:
        # Interpolate. Will return the original geometries for between = 0
        geoms = interpolate_all(geoms, kind=interpol_key, **interpol_kwargs)
        dump_geoms(geoms, "interpolated")

    # Recreate geometries with desired coordinate type and keyword arguments
    geoms = standardize_geoms(geoms, coord_type, geom_kwargs, union=union)

    # Create COS objects and supply a function that yields new Calculators,
    # as needed for growing COS classes, where images are added over time.
    if run_dict["cos"]:
        cos_cls = COS_DICT[cos_key]
        if issubclass(cos_cls, GrowingChainOfStates) or isinstance(
            cos_cls, type(FreezingString)
        ):
            cos_kwargs["calc_getter"] = get_calc_closure("image", calc_key, calc_kwargs)
        geom = COS_DICT[cos_key](geoms, **cos_kwargs)
    elif len(geoms) == 1:
        geom = geoms[0]

    if run_dict["stocastic"]:
        stoc_kwargs["calc_kwargs"] = calc_kwargs
        stocastic = STOCASTIC_DICT[stoc_key](geom, **stoc_kwargs)
        stocastic = run_stocastic(stocastic)
    elif run_dict["md"]:
        md_kwargs = run_dict["md"].copy()
        run_md(geom, calc_getter, md_kwargs)
    elif run_dict["scan"]:
        scan_kwargs = run_dict["scan"]
        scan_geoms, scan_vals, scan_energies = run_scan(geom, calc_getter, scan_kwargs)
    elif run_dict["perf"]:
        perf_results = run_perf(geom, calc_getter, **run_dict["perf"])
        print_perf_results(perf_results)
    # This case will handle most pysisyphus runs. A full run encompasses
    # the following steps:
    #
    #    (0. Preoptimization, already handled)
    #     1. (COS)-Optimization
    #     2. TS-Optimization by TSHessianOptimizer or Dimer method
    #     3. IRC integration
    #     4. Optimization of IRC endpoints
    #
    # Everything can be chained. All functions operate on the 'geom' object,
    # which is propagated along through all functions calls.
    #
    # All keys are present in 'run_dict', but most of the corresponding values will
    # be set to zero.
    elif any(
        [run_dict[key] is not None for key in ("opt", "tsopt", "irc", "mdp", "endopt")]
    ):

        #######
        # OPT #
        #######

        if run_dict["opt"]:
            if run_dict["shake"]:
                print(highlight_text("Shake coordinates"))
                shaked_coords = shake_coords(geom.coords, **run_dict["shake"])
                geom.coords = shaked_coords
                print(f"Shaken coordinates:\n{geom.as_xyz()}")
            opt_result = run_opt(
                geom, calc_getter, opt_key, opt_kwargs, print_thermo=True
            )
            opt_geom = opt_result.geom
            opt = opt_result.opt
            # Keep a backup of the optimized geometry
            if isinstance(opt_geom, ChainOfStates.ChainOfStates):
                # Set some variables that are later collected into RunResult
                cos = opt_geom
                cos_opt = opt
                # copy() is not present for ChainOfState objects, so we just keep
                # using the COS object with a different name.
                geom = opt_geom
            else:
                geom = opt_geom.copy()

        #########
        # TSOPT #
        #########

        abort = False
        if run_dict["tsopt"]:
            # Use a separate implementation for TS-Optimizations started from
            # COS-optimizations.
            try:
                if isinstance(geom, ChainOfStates.ChainOfStates):
                    ts_calc_getter = get_calc_closure(tsopt_key, calc_key, calc_kwargs)
                    ts_opt_result = run_tsopt_from_cos(
                        geom, tsopt_key, tsopt_kwargs, ts_calc_getter
                    )
                else:
                    ts_opt_result = run_opt(
                        geom,
                        calc_getter,
                        tsopt_key,
                        tsopt_kwargs,
                        title="TS-Optimization",
                        copy_final_geom="ts_opt.xyz",
                    )
                ts_geom = ts_opt_result.geom
                ts_opt = ts_opt_result.opt
                geom = ts_geom.copy_all()
            except HEIIsFirstOrLastException:
                abort = True

        #######
        # IRC #
        #######

        ran_irc = False
        if (not abort) and run_dict["irc"]:
            # After a Dimer run we continue with the actual calculator
            # and not the Dimer calculator.
            if calc_key == "dimer":
                calc_getter = act_calc_getter
            irc_geom = geom.copy()
            irc = run_irc(geom, irc_key, irc_kwargs, calc_getter)
            # IRC geom won't have a calculator, so we set the appropriate values here.
            irc_geom.energy = irc.ts_energy
            irc_geom.cart_hessian = irc.init_hessian
            ran_irc = True

        ##########
        # ENDOPT #
        ##########

        # Run 'endopt' when a previous IRC calculation was done
        # if ran_irc and run_dict["endopt"]:
        if run_dict["endopt"]:
            if not ran_irc:

                _, irc_geom, _ = geoms  # IRC geom should correspond to the TS

                class DummyIRC:
                    def __init__(self, geoms):
                        first, _, last = geoms
                        self.atoms = copy.copy(first.atoms)
                        self.forward = True
                        self.backward = True
                        self.downhill = False
                        self.all_coords = [
                            geom.cart_coords.copy() for geom in (first, last)
                        ]
                        self.hessian_init = "dummy"

                irc = DummyIRC(geoms)

            do_thermo = run_dict["endopt"].get("do_hess", False) and can_thermoanalysis
            T = run_dict["endopt"]["T"]
            p = run_dict["endopt"]["p"]
            # Order is forward, backward, downhill
            endopt_results = run_endopt(irc, endopt_key, endopt_kwargs, calc_getter)

            # Determine "left" and "right" geoms
            # Only downhill
            if irc.downhill:
                left_results = endopt_results[2]
            # Only backward
            elif irc.backward and not irc.forward:
                left_results = endopt_results[1]
            # Forward and backward run
            else:
                left_results = endopt_results[0]

            left_geoms = [result.geom for result in left_results]
            left_fns = [result.fn for result in left_results]
            # Use 'backward' results; might be empty.
            right_geoms = [result.geom for result in endopt_results[1]]
            right_fns = [result.fn for result in endopt_results[1]]

            end_geoms = left_geoms + right_geoms
            if run_dict["cos"]:
                end_fns = left_fns + right_fns
                do_rmsds(xyz, geoms, end_fns, end_geoms)

            # Try to compute barriers.
            barrier_kwargs = {
                "do_thermo": do_thermo,
                "T": T,
                "p": p,
                "calc_getter": calc_getter,
                "solv_calc_getter": solv_calc_getter,
            }
            barrier_kwargs_ = run_dict.get("barriers", {})
            barrier_kwargs.update(barrier_kwargs_)
            do_endopt_ts_barriers(
                irc_geom,
                left_geoms,
                right_geoms,
                left_fns=left_fns,
                right_fns=right_fns,
                **barrier_kwargs,
            )

            # Dump TS and endopt geoms to .trj. But only when we did not optimize
            # separate fragments.
            if len(left_geoms) == 1 and len(right_geoms) in (0, 1):
                trj_fn = "left_ts_right_geoms.trj"
                write_geoms_to_trj(
                    list(it.chain(left_geoms, [irc_geom], right_geoms)),
                    trj_fn,
                    comments=list(it.chain(left_fns, ["IRC start"], right_fns)),
                )
                print(f"Wrote optimized end-geometries and TS to '{trj_fn}'")

        if run_dict["mdp"]:
            mdp_kwargs = run_dict["mdp"]
            mdp_result = run_mdp(geom, calc_getter, mdp_kwargs)
    # Fallback when no specific job type was specified
    else:
        calced_geoms, calced_results = run_calculations(
            geoms, calc_getter, scheduler, run_func=calc_run_func
        )

    # We can't use locals() in the dict comprehension, as it runs in its own
    # local scope.
    locals_ = locals()
    results = {key: locals_.get(key, None) for key in RunResult._fields}
    run_result = RunResult(**results)
    return run_result


def check_asserts(results, run_dict):
    print(highlight_text(f"Asserting results"))

    assert_ = run_dict["assert"]
    keys = list(assert_.keys())
    objs_attrs = [key.split(".") for key in keys]
    ref_vals = [assert_[k] for k in keys]
    matches = list()
    for i, ((obj, attr), ref_val) in enumerate(zip(objs_attrs, ref_vals)):
        cur_val = getattr(getattr(results, obj), attr)
        matched = approx_float(cur_val, ref_val)
        print(f"{i:02d}: {obj}.{attr}")
        print(f"\tReference: {ref_val}")
        print(f"\t  Current: {cur_val}")
        print(f"\t  Matches: {bool_color(matched)}")
        matches.append(matched)

    assert all(matches)
    print()


def do_clean(force=False):
    """Deletes files from previous runs in the cwd.
    A similar function could be used to store everything ..."""
    cwd = Path(".").resolve()
    rm_globs = (
        "cycle*.trj",
        "interpolated.trj",
        "interpolated.image*.xyz",
        "calculator.log",
        "optimizer.log",
        "tsoptimizer.log",
        "wfoverlap.log",
        "host_*.calculator.log",
        "host_*.wfoverlap.log",
        "wfo_*.out" "optimization.trj",
        "cos.log",
        "*.gradient",
        "optimizer_results.yaml",
        # ORCA specific
        "*.orca.gbw",
        "*.orca.cis",
        "*.orca.engrad",
        "*.orca.hessian",
        "*.orca.inp",
        "*.orca.hess",
        "*.orca.molden",
        # OpenMOLCAS specific
        "calculator*.out",
        "calculator*.JobIph",
        "calculator*.RasOrb",
        "*rasscf.molden",
        # Turbomole specific
        "calculator_*.control",
        "calculator_*.coord",
        "calculator_*.mos",
        "calculator_*.ciss_a",
        "calculator*.sing_a",
        "*wavefunction.molden",
        "*input.xyz",
        "*.coord",
        # PySCF specific
        "calculator*.chkfile",
        "*.pyscf.out",
        "*.chkfile",
        # WFOverlap specific
        "wfo_*.*.out",
        # XTB specific
        "image*.grad",
        "calculator*.grad",
        "calculator*.xcontrol",
        "calculator*.charges",
        "image_*",
        "splined_ts_guess.xyz",
        "splined_hei_tangent",
        "cart_hei_tangent.trj",
        "dimer_ts.xyz",
        "dimer_pickle",
        "interpolated.geom_*.xyz",
        # Wavefunction overlap
        "wfo_*",
        "image*.molden",
        "jmol.spt",
        "overlap_data.h5",
        "*_CDD.png",
        "*_CDD.cub",
        "internal_coords.log",
        "hei_tangent",
        "optimization.trj",
        "splined_hei.xyz",
        "ts_opt.xyz",
        "final_geometry.xyz",
        "calculated_init_hessian",
        "cur_out",
        # HDF5 files
        "optimization.h5",
        "afir.h5",
        # Optimization files
        "*_optimization.trj",
        # Preopt files
        "first_*",
        "last_*",
        # TSOpt
        "rsirfo*",
        # IRC files
        "irc_*",
        "irc.log",
        "finished_*",
        # IRC/Endopt files
        "backward_*",
        "forward_*",
        # Misc
        "*imaginary_mode_*.trj",
        "cart_hei_tangent",
        "ts_calculated_init_cart_hessian",
        "calculated_final_cart_hessian",
        "*final_geometry.xyz",
        "*final_geometries.trj",
        "current_geometry.xyz",
        "*current_geometries.trj",
        "hess_calc_cyc*.h5",
        "ts_hess_calc_cyc*.h5",
        "hess_init_irc.h5",
        "final_hessian.h5",
        "ts_current_geometry.xyz",
        "dimer_*",
        "plain_hei_tangent",
        "plain_hei.xyz",
        "hess_calc_irc*.h5",
        "rebuilt_primitives.xyz",
        "RUN.yaml",
        "middle_for_preopt.trj",
        "relaxed_scan.trj",
        "too_similar.trj",
        # MDP
        "mdp_ee_ascent.trj",
        "mdp_ee_fin_*.trj",
        "mdp_ee_init_*.trj",
        "aligned.geom*xyz",
        "cos_hei.trj",
        # Dimer
        "calculator_*.N",
        "calculator_*.N.trj",
        "dimer.log",
        "*.gfnff_topo",
        # DFTB+
        "*.detailed.out",
        "*.geometry.gen",
        "*.dftb_in.hsd",
        "*.EXC.DAT",
        "*.XplusY.DAT",
        "*.dftb.out",
        "rsprfo_*",
        "reparametrized.trj",
        "end_geoms_and_ts.trj",
        "left_ts_right_geoms.trj",
        "ts_final_hessian.h5",
        "third_deriv.h5",
        "*.ao_ovlp_rec",
        # MOPAC
        "*.mopac.aux",
        "*.mopac.arc",
        "*.mopac.mop",
        "*.mopac.out",
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
                print(f"Deleted {p}")
            except FileNotFoundError:
                pass
        try:
            os.unlink("cur_out")
        except FileNotFoundError:
            pass

    if force:
        delete()
        return
    # If we dont force the cleaning ask for confirmation first
    elif to_rm_paths and confirm_input("Delete these files?"):
        delete()
    else:
        print("No files found for removal.")


def print_header():
    """Generated from https://asciiartgen.now.sh/?s=pysisyphus&style=colossal"""
    normal_logo = """                           d8b                            888
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

    xmas_logo = r"""
                                   \|/
                           d8b    --o--                   888              X
                  x        Y8P     /|\                    888
                                                          888
88888b.  888  888 .d8888b  888 .d8888b  888  888 88888b.  88888b.  888  888 .d8888b
888 "88b 888  888 88K      888 88K      888  888 888 "88b 888 "88b 888  888 88K
888  888 888  888 "Y8888b. 888 "Y8888b. 888  888 888  888 888  888 888  888 "Y8888b.
888 d88P Y88b 888      X88 888      X88 Y88b 888 888 d88P 888  888 Y88b 888      X88
88888P"   "Y88888  88888P' 888  88888P'  "Y88888 88888P"  888  888  "Y88888  88888P'
888  |        888   |  |    |     |          888 888       |    |              |
888  O   Y8b d88P   o  |    X    / \    Y8b d88P 888       x   / \             O
888       "Y88P"       O         \_/     "Y88P"  888           \_/              
            |                         x           |                      \|/
   X        |                        xox          O                     --X--
            O                         x                                  /|\\
"""
    now = datetime.datetime.now()
    today = now.date()
    xmas = (today.month == 12) and (today.day <= 24)
    logo = xmas_logo if xmas else normal_logo
    version = f"Version {__version__}"
    try:
        commit = re.compile(r"g(\w+)\.").search(version).group(1)
        commit_line = f"Git commit {commit}\n"
    # Raised when regex search failed
    except AttributeError:
        commit_line = ""
    vi = sys.version_info
    sv = f"{vi.major}.{vi.minor}.{vi.micro}"  # Python
    npv = np.__version__  # Numpy
    spv = sp.__version__  # SciPy
    cwd = Path(".").resolve()
    print(
        f"{logo}\n\n{version} (Python {sv}, NumPy {npv}, SciPy {spv})\n"
        f"{commit_line}"
        f"Executed at {now.strftime('%c')} on '{platform.node()}'\n"
        f"Platform: {platform.platform()}\n"
        f"Interpreter: {sys.executable}\n"
        f"Current working directory: {cwd}\n"
    )


def print_bibtex():
    bibtex = textwrap.dedent(
        """@article{Steinmetzer2020,
        doi = {10.1002/qua.26390},
        url = {https://doi.org/10.1002/qua.26390},
        year = {2020},
        month = aug,
        publisher = {Wiley},
        author = {Johannes Steinmetzer and Stephan Kupfer and Stefanie Gräfe},
        title = {pysisyphus: Exploring potential energy surfaces in ground and excited states},
        journal = {International Journal of Quantum Chemistry}
    }"""
    )
    print(bibtex)


def run_from_dict(
    run_dict,
    cwd=None,
    set_defaults=True,
    yaml_fn=None,
    cp=None,
    scheduler=None,
    clean=False,
    fclean=False,
    version=False,
    restart=False,
):
    if cwd is None:
        cwd = Path(".")

    print_header()

    # Citation
    citation = (
        "If pysisyphus benefitted your research please cite:\n\n"
        "\thttps://doi.org/10.1002/qua.26390\n\nGood luck!\n"
    )
    print(citation)

    init_logging(cwd, scheduler)
    # Load defaults etc.
    if set_defaults:
        run_dict = setup_run_dict(run_dict)
    sys.stdout.flush()

    if cp:
        copy_yaml_and_geometries(run_dict, yaml_fn, cp)
        return
    elif clean:
        do_clean()
        return
    elif fclean:
        do_clean(force=True)
        return
    # Return after header was printed
    elif version:
        return

    run_dict_without_none = {k: v for k, v in run_dict.items() if v is not None}
    pprint(run_dict_without_none)
    print()
    sys.stdout.flush()

    run_result = main(run_dict, restart, cwd, scheduler)

    if run_dict["assert"] is not None:
        print()
        check_asserts(run_result, run_dict)

    return run_result


def run():
    start_time = datetime.datetime.now()
    args = parse_args(sys.argv[1:])

    # Defaults
    run_dict = {}
    yaml_dir = Path(".")

    if args.yaml:
        yaml_dir = Path(os.path.abspath(args.yaml)).parent
        with open(args.yaml) as handle:
            yaml_str = handle.read()
        try:
            loader = get_loader()
            try:
                run_dict = yaml.load(yaml_str, Loader=loader)
            except yaml.constructor.ConstructorError as err:
                mobj = re.compile("for the tag '\!(\w+)'").search(err.problem)
                if mobj:
                    err_unit = mobj.group(1)
                    best_match, _ = find_closest_sequence(err_unit, UNITS)
                    print(
                        f"Unknown unit!\nKnown units are\n'{UNITS}'.\n"
                        f"Did you mean '{best_match}', instead of '{err_unit}'?\n"
                    )
                raise err
            assert type(run_dict) == type(dict())
        except (AssertionError, yaml.parser.ParserError) as err:
            print(err)
            if not (args.yaml.lower().endswith(".yaml")):
                print("Are you sure that you supplied a YAML file?")
            sys.exit(1)
    elif args.bibtex:
        print_bibtex()
        return

    run_kwargs = {
        "cwd": yaml_dir,
        "set_defaults": True,
        "yaml_fn": args.yaml,
        "cp": args.cp,
        "scheduler": args.scheduler,
        "clean": args.clean,
        "fclean": args.fclean,
        "version": args.version,
        "restart": args.restart,
    }
    run_result = run_from_dict(run_dict, **run_kwargs)

    end_time = datetime.datetime.now()
    duration = end_time - start_time
    # Only keep hh:mm:ss
    duration_hms = str(duration).split(".")[0]
    print(f"pysisyphus run took {duration_hms} h.")

    return 0


if __name__ == "__main__":
    run()
