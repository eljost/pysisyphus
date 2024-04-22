# [1] https://doi.org/10.1021/acs.jctc.0c01306
#     Single-Point Hessian Calculations for Improved Vibrational Frequencies and
#     Rigid-Rotor-Harmonic-Oscillator Thermodynamics
#     Spicher, Grimme


from dataclasses import dataclass
from math import floor, ceil
from pathlib import Path
import shutil
from typing import Callable, Dict, Optional, Tuple
import sys

import numpy as np

from pysisyphus.calculators import Dimer, ExternalPotential
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
    # Converge the lowest two modes for TS optimizations, as only one must have a
    # negative eigenvalue. For minimizations the lowest mode is enough.
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

        # Initialize dask cluster, if required
        cluster = geom.init_dask()

        try:
            if geom.images[0].calculator.track:
                print(
                    "Propagating root information of first image along COS w/ "
                    "energy calculations."
                )
                geom.propagate()
        except AttributeError:
            pass
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
        print(highlight_text(f"Running {title}", level=level) + "\n")
        print(f"     Input geometry: {geom.describe()}")
        print(f"  Coordinate system: {geom.coord_type}")
        print(f"  Coordinate number: {len(geom.coords)}")
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
                print("\tcalculated wavefunction for image {j:03d}")
                cur_calc = image.calculator
                try:
                    next_calc = geom.images[j + 1].calculator
                except IndexError:
                    pass
                try:
                    next_calc.set_chkfiles(cur_calc.get_chkfiles())
                    print("\tset wavefunction of image {j:03d} on image {j+1:03d}")
                except AttributeError:
                    break
                sys.stdout.flush()

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

    if is_cos:
        # Shutdown dask cluster, if required
        geom.exit_dask(cluster)

    opt_result = OptResult(opt, opt.geometry, opt.final_fn)
    return opt_result


def get_optimal_bias(
    ref_geom: Geometry,
    calc_getter: Callable,
    opt_key: str,
    opt_kwargs: Dict,
    k_max: float,
    k_min: float = 0.0,
    rmsd_target: float = 0.188973,
    rmsd_thresh: Optional[float] = None,
    rmsd_kwargs: Optional[Dict] = None,
    k_thresh: float = 1e-3,
    strict=True,
) -> Tuple[OptResult, float, bool]:
    """Driver to determine optimal bias value k for RMSD restraint,
    as required in single point hessian (SPH) calculations.

    Parameters
    ----------
    ref_geom
        Reference geometry. Starting point of the optimizations and reference
        for RMSD calculations.
    calc_getter
        Function that returns the actual calculator, providing the energy and
        its derivatives.
    opt_key
        Determines optimizer type. See pysisyphus.optimizers.cls_map.
    opt_kwargs
        Optional dict of arguments passed to the optimizer.
    k_max
        Maximum absolute value of bias factor k. Must be a > k_min.
    k_min
        Minimum absolute value of bias factor k. Must be a positive number >= 0.0.
        Defaults to 0.0.
    rmsd_target
        Target RMSD value in au. Defaults to 0.188973 a0 (approx. 0.1 Å).
    rmsd_thresh
        Allowed deviation from rmsd_target in au. If omitted, 5% of rmsd_target are used.
    rmsd_kwargs
        Additional keyword arguments that are passed to the RMSD class, e.g., atom_indices.
    k_thresh:
       When the absolute value of k_bias - k_min or k_max becomes smaller than
       k_thresh, the bisection is aborted.
    strict
        If True, AssertionError is raised when an optimization did not converged.

    Returns
    -------
    opt_result
        OptimizationResult object containig the Optimizer object.
    k_opt
        Optimal value of k_bias.
    valid_k
        Whether an appropriate bias value k was found.
    """
    assert rmsd_target > 0.0
    assert k_min >= 0.0
    assert k_max > k_min
    assert k_thresh >= 0.0
    if rmsd_thresh is None:
        rmsd_thresh = 0.05 * rmsd_target
    if rmsd_kwargs is None:
        rmsd_kwargs = dict()

    opt_counter = 0

    def run_biased_opt(k):
        nonlocal opt_counter
        geom = ref_geom.copy()

        def biased_calc_getter():
            act_calc = calc_getter()
            _rmsd_kwargs = {
                "type": "rmsd",
                "k": k,
            }
            _rmsd_kwargs.update(rmsd_kwargs)
            potentials = (_rmsd_kwargs,)
            calc = ExternalPotential(
                calculator=act_calc, potentials=potentials, geom=ref_geom
            )
            return calc

        sys.stdout.flush()
        prefix = f"biased_{opt_counter:02d}"
        prefixed_opt_kwargs = opt_kwargs.copy()
        prefixed_opt_kwargs.update(
            {
                "prefix": prefix,
                "h5_group_name": f"{prefix}_opt",
            }
        )
        print(f"@@@ Running optimization {opt_counter:02d} with {k=:.6f} au.")
        opt_result = run_opt(geom, biased_calc_getter, opt_key, prefixed_opt_kwargs)
        if strict:
            assert opt_result.opt.is_converged

        rmsd_current = ref_geom.rmsd(geom)
        print(
            f"@@@ Biased optimization {opt_counter:02d} with {k=:.6f} au yielded an "
            f"RMSD of RMSD={rmsd_current:.4f} au."
        )
        opt_counter += 1
        return opt_result, rmsd_current

    k_bias = 0.0
    opt_result, rmsd_current = run_biased_opt(
        k=k_bias
    )  # Start with unbiased optimization
    unbiased_failed = rmsd_current > rmsd_target

    k_min0 = k_min
    k_max0 = k_max

    def k_is_valid(k):
        """Return whether k is too close to either k_min0 or k_max0."""
        return all([abs(k_bias - k) >= k_thresh for k in (k_min0, k_max0)])

    valid_k = True

    # Only run the loop when the initial unbiased optimization failed. This way
    # we don't have to check for an early return before the loop, but can keep
    # one return statement at the end of the function.
    while unbiased_failed and abs(rmsd_current - rmsd_target) > rmsd_thresh:
        k_bias = -0.5 * (k_min + k_max)
        if not (valid_k := k_is_valid(k_bias)):
            break
        # Biased geometry optimization
        opt_result, rmsd_current = run_biased_opt(k_bias)
        if rmsd_current < rmsd_target:
            k_max = abs(k_bias)
        elif rmsd_current > rmsd_target:
            k_min = abs(k_bias)
        print(
            f"@@@ After biased optimization {opt_counter:02d}: {k_min=:.4f}, {k_bias=:.4f}, "
            f"{k_max=:.4f}, {rmsd_current=:4f} au."
        )
    k_opt = k_bias
    return opt_result, k_opt, valid_k
