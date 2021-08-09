from math import floor, ceil
import shutil

import numpy as np

from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import do_final_hessian
from pysisyphus.helpers_pure import highlight_text, report_frozen_atoms
from pysisyphus.modefollow import NormalMode, geom_davidson
from pysisyphus.optimizers import *
from pysisyphus.optimizers.hessian_updates import bfgs_update
from pysisyphus.tsoptimizers import *

OPT_DICT = {
    "bfgs": BFGS.BFGS,
    "cg": ConjugateGradient.ConjugateGradient,
    "fire": FIRE.FIRE,
    "lbfgs": LBFGS.LBFGS,
    "micro": MicroOptimizer,
    "nc": NCOptimizer.NCOptimizer,
    "oniom": ONIOMOpt,
    "plbfgs": PreconLBFGS.PreconLBFGS,
    "psd": PreconSteepestDescent.PreconSteepestDescent,
    "qm": QuickMin.QuickMin,
    "rfo": RFOptimizer.RFOptimizer,
    "sd": SteepestDescent.SteepestDescent,
    "sqnm": StabilizedQNMethod.StabilizedQNMethod,
    "string": StringOptimizer.StringOptimizer,
}

TSOPT_DICT = {
    "dimer": dimer_method,
    "rsprfo": RSPRFOptimizer,
    "trim": TRIM,
    "rsirfo": RSIRFOptimizer,
}


def get_opt_cls(opt_key):
    try:
        opt_cls = OPT_DICT[opt_key]
    except KeyError:
        opt_cls = TSOPT_DICT[opt_key]
    return opt_cls


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




def run_opt(
    geom,
    calc_getter,
    opt_key,
    opt_kwargs,
    title="Optimization",
    copy_final_geom=None,
    level=0,
):
    is_cos = issubclass(type(geom), ChainOfStates)

    if is_cos:
        for image in geom.images:
            image.set_calculator(calc_getter())
            title = str(geom)
    else:
        geom.set_calculator(calc_getter())

    do_hess = opt_kwargs.pop("do_hess", False)
    do_davidson = opt_kwargs.pop("do_davidson", False)
    T = opt_kwargs.pop("T", 298.15)

    opt = get_opt_cls(opt_key)(geom, **opt_kwargs)
    print(highlight_text(f"Running {title}", level=level))
    print(f"\nInput structure: {geom.describe()}")
    report_frozen_atoms(geom)
    print()

    opt.run()

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
        do_final_hessian(geom, write_imag_modes=True, prefix=prefix, T=T)
    print()

    return opt.geometry, opt
