from dataclasses import dataclass
from math import floor, ceil
from pathlib import Path
import shutil
import sys

import numpy as np

from pysisyphus.calculators import Dimer
from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.config import T_DEFAULT, p_DEFAULT
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import do_final_hessian
from pysisyphus.helpers_pure import highlight_text, report_frozen_atoms
from pysisyphus.io import save_hessian
from pysisyphus.modefollow import NormalMode, geom_davidson
from pysisyphus.optimizers.cls_map import get_opt_cls, key_is_tsopt
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer
from pysisyphus.optimizers.hessian_updates import bfgs_update
from pysisyphus.tsoptimizers.TSHessianOptimizer import TSHessianOptimizer


def opt_davidson(opt, tsopt=True, res_rms_thresh=1e-4):
    try:
        H = opt.H
    except AttributeError:
        if tsopt:
            raise Exception("Can't handle TS optimization without Hessian yet!")

        # Create approximate updated Hessian
        cart_coords = opt.cart_coords
        cart_forces = opt.cart_forces
        coord_diffs = np.diff(cart_coords, axis=0)
        grad_diffs = -np.diff(cart_forces, axis=0)
        H = np.eye(cart_coords[0].size)
        for s, y in zip(coord_diffs, grad_diffs):
            dH, _ = bfgs_update(H, s, y)
            H += dH

    geom = opt.geometry
    if geom.coord_type != "cart":
        H = geom.internal.backtransform_hessian(H)

    masses_rep = geom.masses_rep
    # Mass-weigh and project Hessian
    H = geom.eckart_projection(geom.mass_weigh_hessian(H))
    w, v = np.linalg.eigh(H)
    inds = [0, 1] if tsopt else [6]
    # Converge the lowest two modes for TS optimizations; for minimizations the lowest
    # mode would is enough.
    lowest = 2 if tsopt else 1
    guess_modes = [NormalMode(l, masses_rep) for l in v[:, inds].T]

    davidson_kwargs = {
        "hessian_precon": H,
        "guess_modes": guess_modes,
        "lowest": lowest,
        "res_rms_thresh": res_rms_thresh,
        "remove_trans_rot": True,
    }

    result = geom_davidson(geom, **davidson_kwargs)
    return result


@dataclass
class OptResult:
    opt: Optimizer
    geom: Geometry
    fn: Path


def run_opt(
    geom,
    calc_getter,
    opt_key,
    opt_kwargs=None,
    iterative=False,
    iterative_max_cycles=5,
    iterative_thresh=-15,
    iterative_scale=2.00,
    cart_hessian=None,
    print_thermo=False,
    title="Optimization",
    copy_final_geom=None,
    level=0,
):
    is_cos = issubclass(type(geom), ChainOfStates)
    is_tsopt = key_is_tsopt(opt_key)
    # Disallow iterative optimizations for COS objects
    if opt_kwargs is None:
        opt_kwargs = dict()
    is_iterative = (not is_cos) and (iterative or opt_kwargs.pop("iterative", False))

    if is_cos:
        # Set calculators on all images
        for image in geom.images:
            image.set_calculator(calc_getter())
            title = str(geom)
    else:
        geom.set_calculator(calc_getter())
        geom.cart_hessian = cart_hessian

    do_hess = opt_kwargs.pop("do_hess", False)
    do_davidson = opt_kwargs.pop("do_davidson", False)
    T = opt_kwargs.pop("T", T_DEFAULT)
    p = opt_kwargs.pop("p", p_DEFAULT)
    propagate = opt_kwargs.pop("propagate", False)

    opt_cls = get_opt_cls(opt_key)
    for i in range(iterative_max_cycles):
        # Modify hessian_init in later cycles, te reuse the calculated Hessian
        # from the final geometry of the previous optimization.
        if (i > 0) and issubclass(opt_cls, HessianOptimizer):
            opt_kwargs["hessian_init"] = "calc"
        opt = opt_cls(geom, **opt_kwargs)
        print(highlight_text(f"Running {title}", level=level))
        print(f"\n   Input geometry: {geom.describe()}")
        print(f"  Coordinate system: {geom.coord_type}")
        print(f"         Calculator: {geom.calculator}")
        print(f"             Charge: {geom.calculator.charge}")
        print(f"       Multiplicity: {geom.calculator.mult}")
        print(f"          Optimizer: {opt_key}\n")
        report_frozen_atoms(geom)
        print()

        # Try to propagate chkfiles along calculators in COS optimizations
        if propagate and is_cos:
            print("Propagating chkfiles along COS")
            for j, image in enumerate(geom.images):
                image.energy
                cur_calc = image.calculator
                try:
                    next_calc = geom.images[j + 1].calculator
                except IndexError:
                    pass
                try:
                    next_calc.set_chkfiles(cur_calc.get_chkfiles())
                except AttributeError:
                    break

        opt.run()

        # Only do 1 cycle in non-iterative optimizations
        if not is_iterative:
            break

        # Determine imaginary modes for subsequent displacements
        nus, *_, cart_displs = geom.get_normal_modes()
        below_thresh = nus < iterative_thresh
        # Never displace along transition vector in ts-optimizations. Just skip it.
        if is_tsopt:
            below_thresh[0] = False
        imag_nus = nus[below_thresh]
        imag_displs = cart_displs[:, below_thresh].T

        if len(imag_nus) == 0:
            print(f"Iterative optimization converged in cycle {i}.")
            break

        h5_fn = f"hess_calc_iter_{i:02d}.h5"
        save_hessian(h5_fn, geom)
        print(f"Saved HDF5-Hessian to {h5_fn}.")

        print(f"\nImaginary modes below threshold of {iterative_thresh:.2f} cm⁻¹:")
        for j, nu in enumerate(nus[below_thresh]):
            print(f"\t{j:02d}: {nu:8.2f} cm⁻¹")
        sys.stdout.flush()

        print(f"\nGeometry after optimization cycle {i}:\n{geom.as_xyz()}")

        # Displace along imaginary modes
        for j, (nu, imag_displ) in enumerate(zip(imag_nus, imag_displs)):
            step = iterative_scale * imag_displ
            new_cart_coords = geom.cart_coords + step
            geom.cart_coords = new_cart_coords

        print(f"\nDisplaced geometry for optimization cycle {i+1}:\n{geom.as_xyz()}\n")

    # ChainOfStates specific
    if is_cos and (not opt.stopped):
        hei_coords, hei_energy, hei_tangent, hei_frac_index = geom.get_splined_hei()
        floor_ind = floor(hei_frac_index)
        ceil_ind = ceil(hei_frac_index)
        print(
            f"Splined HEI is at {hei_frac_index:.2f}/{len(geom.images)-1:.2f}, "
            f"between image {floor_ind} and {ceil_ind} (0-based indexing)."
        )
        hei_geom = Geometry(geom.images[0].atoms, hei_coords)
        hei_geom.energy = hei_energy
        hei_fn = "splined_hei.xyz"
        with open(hei_fn, "w") as handle:
            handle.write(hei_geom.as_xyz())
        print(f"Wrote splined HEI to '{hei_fn}'")

    if copy_final_geom and opt.is_converged:
        copy_fn = copy_final_geom
        shutil.copy(opt.final_fn, copy_fn)
        print(f"Copied '{opt.final_fn}' to '{copy_fn}'.")

    if do_davidson and (not opt.stopped):
        tsopt = isinstance(opt, TSHessianOptimizer.TSHessianOptimizer) or isinstance(
            geom.calculator, Dimer
        )
        type_ = "TS" if tsopt else "minimum"
        print(highlight_text(f"Davidson after {type_} search", level=1))
        opt_davidson(opt, tsopt=tsopt)
    elif do_hess and (not opt.stopped):
        print()
        prefix = opt_kwargs.get("prefix", "")
        out_dir = opt_kwargs.get("out_dir", None)
        do_final_hessian(
            geom,
            write_imag_modes=True,
            prefix=prefix,
            T=T,
            p=p,
            print_thermo=print_thermo,
            is_ts=is_tsopt,
            out_dir=out_dir,
        )
    print()

    opt_result = OptResult(opt, opt.geometry, opt.final_fn)
    return opt_result
