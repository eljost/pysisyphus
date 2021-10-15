#!/usr/bin/env python

import argparse
from dataclasses import dataclass
import os
from pathlib import Path
import shutil
import sys
import time
import traceback

import orjson
import numpy as np
import psutil

from pysisyphus.calculators import XTB
from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.cos.GrowingNT import GrowingNT
from pysisyphus.helpers import geom_loader, do_final_hessian, check_for_end_sign
from pysisyphus.helpers_pure import highlight_text, touch
from pysisyphus.irc import EulerPC
from pysisyphus.optimizers.PreconLBFGS import PreconLBFGS
from pysisyphus.tsoptimizers.RSPRFOptimizer import RSPRFOptimizer


def dump_trj(geoms, fn):
    trj = "\n".join([geom.as_xyz() for geom in geoms])
    with open(fn, "w") as handle:
        handle.write(trj)


@dataclass
class CGNTResult:

    atoms: tuple
    init_energy: float
    gnt_last_image_c3d: np.ndarray

    gnt_ts_c3d: np.ndarray = None
    ts_c3d: np.ndarray = None
    ts_opt_converged: bool = False
    final_energy: float = None
    barrier: float = None


def run_gnt(id_, geom, bonds, flag="CONVERGED", opt_ts=False, irc=False, force=False):
    work_dir = Path(id_)

    def get_fn(fn):
        return work_dir / fn

    # Return quick
    flag_fn = get_fn(flag)
    if not force and flag_fn.exists():
        print(f"@@@ {id_} is already converged!")
        return True

    # Else, start over
    try:
        shutil.rmtree(work_dir)
    except FileNotFoundError:
        pass

    os.mkdir(work_dir)

    calc = XTB(
        pal=psutil.cpu_count(logical=False),
        quiet=True,
        gfn=2,
        out_dir=work_dir,
        retry_etemp=1000.0,
    )
    geom.set_calculator(calc)

    # Calculate initial energy and dump initial geometry
    init_energy = geom.energy
    with open(get_fn("input.xyz"), "w") as handle:
        handle.write(geom.as_xyz())

    #
    # Run Growing Newton trajectory
    #
    gnt_kwargs = {
        "step_len": 0.1,
        "bonds": bonds,
        "stop_after_ts": True,
        "rms_thresh": 0.003,
    }
    gnt = GrowingNT(geom, **gnt_kwargs)
    opt_kwargs = {
        "max_cycles": 1_000,
        "dump": True,
        "max_step": 0.1,
        "out_dir": work_dir,
        "line_search": "armijo_fg",
    }
    print(highlight_text(f"{id_} GNT Optimization"))
    opt = PreconLBFGS(gnt, **opt_kwargs)
    try:
        opt.run()
    except Exception:
        exc_info = sys.exc_info()
        traceback.print_exception(*exc_info)
    print()

    # Pick best geometry from GNT. Check if we found a TS (guess). If not
    # use the latest converged geometry instead.
    try:
        last_gnt_geom = gnt.ts_images[0]
    except IndexError:
        last_gnt_geom = gnt.images[-1]
    last_gnt_geom.set_calculator(calc)
    do_final_hessian(
        last_gnt_geom,
        write_imag_modes=True,
        prefix=f"{id_}_last_gnt_geom",
        out_dir=work_dir,
    )
    gnt_succeeded = opt.is_converged

    # Dump all GNT images and all obtained stationary points
    dump_trj(gnt.images, get_fn(f"{id_}_gnt.trj"))
    dump_trj(gnt.sp_images, get_fn(f"{id_}_gnt_sps.trj"))

    # Run actual TS optimization if explicitly requested, or the
    # previous GNT run failed. But don't run TS optimization if
    # the GNT never completed its first step.
    ts_opt = None
    if opt_ts or (not gnt_succeeded and (len(gnt.images) > 1)):
        print(highlight_text(f"{id_} (Recovery) TS optimization"))
        bonds_ = np.array(bonds, dtype=int).reshape(-1, 3)
        typed_prims = [["BOND", from_, to_] for from_, to_, _ in bonds_]

        # Start from latest converged point on Newton trajectory
        ts_guess_geom = last_gnt_geom.copy(
            coord_type="redund",
            # coord_kwargs={
            # "define_prims": typed_prims,
            # },
        )
        ts_guess_geom.set_calculator(calc)
        ts_opt_kwargs = {
            "hessian_recalc": 10,
            "trust_max": 0.5,
            "overachieve_factor": 3,
            "prefix": f"{id_}_ts",
            "out_dir": work_dir,
        }
        ts_opt = RSPRFOptimizer(ts_guess_geom, **ts_opt_kwargs)
        ts_opt.run()
        ts_geom = ts_opt.geometry
        # Overwrite GNT optimizer, so we can just return 'opt.is_converged' as below
        opt = ts_opt
    elif not gnt_succeeded and (len(gnt.images) == 1):
        raise Exception("GNT failed in first step!")
    else:
        ts_geom = last_gnt_geom

    # Below this point we expect a 'ts_geom' to be present
    assert "ts_geom" in locals()

    ts_geom.set_calculator(geom.calculator)
    do_final_hessian(ts_geom, write_imag_modes=True, prefix=id_, out_dir=work_dir)
    final_energy = ts_geom.energy
    barrier = final_energy - init_energy
    barrier_j = barrier * AU2KJPERMOL
    print(f"@@@ {id_}: barrier={barrier_j:.2f} kJ mol⁻¹")

    if opt.is_converged:
        # Touch flag
        touch(flag_fn)

    try:
        gnt_ts_c3d = gnt.ts_images[-1].coords3d
    except IndexError:
        gnt_ts_c3d = None

    cgntres_kwargs = {
        "atoms": geom.atoms,
        "init_energy": init_energy,
        "gnt_last_image_c3d": gnt.images[-1].coords3d,
        "gnt_ts_c3d": gnt_ts_c3d,
        "final_energy": final_energy,
        "barrier": barrier,
    }
    if ts_opt is not None:
        cgntres_kwargs.update(
            ts_c3d=ts_geom.coords3d,
            ts_opt_converged=ts_opt.is_converged,
        )

    if ts_opt and ts_opt.is_converged and irc:
        irc_geom = ts_geom.copy(coord_type="cart")
        irc_geom.set_calculator(calc)
        irc = EulerPC(
            irc_geom,
            out_dir=work_dir,
            hessian_recalc=10,
            prefix=id_,
            hard_rms_grad_thresh=5e-4,
        )
        irc.run()
        print(f"@@@ {id_}: IRC converged? {irc.converged}")

    cgntres = CGNTResult(**cgntres_kwargs)
    with open(get_fn(f"{id_}_dump.json"), "wb") as handle:
        handle.write(
            orjson.dumps(
                cgntres,
                option=orjson.OPT_SERIALIZE_DATACLASS | orjson.OPT_SERIALIZE_NUMPY,
            )
        )

    return opt.is_converged


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("name")
    parser.add_argument("geoms")
    parser.add_argument("bonds", nargs="+", type=int)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--optts", dest="opt_ts", action="store_true")
    parser.add_argument("--irc", action="store_true")
    parser.add_argument("--first", default=None, type=int)
    return parser.parse_args(args)


def run():
    start = time.time()
    args = parse_args(sys.argv[1:])

    geoms = geom_loader(args.geoms)
    print(f"@@@ Loaded {len(geoms)} geometries from '{args.geoms}'")
    if first := args.first:
        geoms = geoms[:first]
        print(f"@@@ Using only first {first} geometries.")
    bonds = np.array(args.bonds, dtype=int).reshape(-1, 3).tolist()
    conv_num = 0
    for i, geom in enumerate(geoms):
        id_ = f"{i:04d}_{args.name}"
        try:
            converged = run_gnt(
                id_, geom, bonds, opt_ts=args.opt_ts, force=args.force, irc=args.irc
            )
        except Exception:
            exc_info = sys.exc_info()
            traceback.print_exception(*exc_info)
            converged = "error"
        print(f"@@@ {id_}: converged? {converged}")
        print("@@@")
        conv_num += 1 if (converged is True) else 0
        if check_for_end_sign():
            break
        print()
    print(f"@@@ converged: {conv_num}/{len(geoms)}")
    dur = time.time() - start
    print(f"@@@ crestnt.py run took {dur/60:.2f} min")


if __name__ == "__main__":
    run()
