#!/usr/bin/env python3

import logging
import os
from pprint import pprint
import sys

import numpy as np

from pysisyphus.helpers import geom_from_library, geom_from_xyz_file
from pysisyphus.calculators.XTB import XTB
from pysisyphus.calculators.Gaussian16 import Gaussian16
from pysisyphus.optimizers.hessian_updates import bfgs_update, flowchart_update
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.RSRFOptimizer import RSRFOptimizer


np.set_printoptions(suppress=True, precision=4)


def rfo(gradient, H, trust=0.3):
    H_aug = np.bmat(
                ((H, gradient[:,None]),
                 (gradient[None,:], [[0]]))
    )
    eigvals, eigvecs = np.linalg.eigh((H_aug+H_aug.T)/2)
    aug_step = np.asarray(eigvecs[:,0]).flatten()
    lambda_ = aug_step[-1]
    step = aug_step[:-1] / lambda_
    norm = np.linalg.norm(step)
    if np.linalg.norm(step) > trust:
        step = step / norm * trust
    return step


def run_bare_rfo(xyz_fn, charge, mult, trust=0.3, max_cycles=150):
    geom = geom_from_library(xyz_fn, coord_type="redund")
    # geom = geom_from_library(xyz_fn)
    geom.set_calculator(XTB(pal=4, charge=charge, mult=mult))
    # geom.set_calculator(Gaussian16(pal=3, route="HF 3-21G", charge=charge, mult=mult))
    grads = list()
    steps = list()
    H = geom.get_initial_hessian()
    converged = False
    for i in range(max_cycles):
        grad = -geom.forces
        grads.append(grad)

        if i > 0:
            dx = steps[-1]
            dg = grads[-1] - grads[-2]
            dH, _ = bfgs_update(H, dx, dg)
            H += dH
        H_proj = H.copy()
        if geom.internal:
            H_proj = geom.internal.project_hessian(H_proj)
        step = rfo(grad, H_proj, trust=trust)
        steps.append(step)
        max_g = np.abs(grad).max()
        rms_g = np.sqrt(np.mean(grad**2))
        max_s = np.abs(step).max()
        rms_s = np.sqrt(np.mean(step**2))
        print(f"{i:02d}: max(f)={max_g:.6f}, rms(f)={rms_g:.6f}, "
              f"max(step)={max_s:.6f}, rms(step)={rms_s:.6f}")
        converged = ((max_g < 4.5e-4) and (rms_g < 1.5e-4)
                     and (max_s < 1.8e-3) and (rms_s < 1.2e-3))
        if converged:
            print("Converged")
            break
        new_coords = geom.coords + step
        geom.coords = new_coords
    return converged, i+1


def run():
    GEOMS = {
        "artemisin.xyz": (0, 1),
        "avobenzone.xyz": (0, 1),
        "azadirachtin.xyz": (0, 1),
        "bisphenol_a.xyz": (0, 1),
        "cetirizine.xyz": (0, 1),
        "codeine.xyz": (0, 1),
        "diisobutyl_phthalate.xyz": (0, 1),
        "easc.xyz": (0, 1),
        "estradiol.xyz": (0, 1),
        "inosine.xyz": (1, 1),
        "maltose.xyz": (0, 1),
        "mg_porphin.xyz": (0, 1),
        "ochratoxin_a.xyz": (0, 1),
        "penicillin_v.xyz": (0, 1),
        "raffinose.xyz": (0, 1),
        "sphingomyelin.xyz": (0, 1),
        "tamoxifen.xyz": (0, 1),
        "vitamin_c.xyz": (0, 1),
        "zn_edta.xyz": (-2, 1)
    }

    results = dict()
    tot_cycles = 0
    fails = 0
    for xyz_fn, (charge, mult) in GEOMS.items():
        full_xyz = "birkholz/" + xyz_fn
        print(f"Running '{full_xyz}'")

        # geom = geom_from_library(full_xyz, coord_type="redund")
        # geom.set_calculator(XTB(pal=2))
        # opt = RFOptimizer(geom, trust_radius=0.5, update_trust=False)
        # opt.run()

        # geom = geom_from_library(full_xyz, coord_type="redund")
        # geom.set_calculator(XTB(pal=2))
        # opt = RSRFOptimizer(geom, trust_radius=0.5, update_trust=True,
                            # hess_update="flowchart", max_micro_cycles=25)
        # opt.run()

        try:
            converged, cycles = run_bare_rfo(full_xyz, charge, mult, trust=0.5)
        except:
            logging.exception("Something went wrong")
            converged = False
            cycles = 0
        # geom = geom_from_library(full_xyz, coord_type="redund")
        # geom.set_calculator(XTB(pal=2, charge=charge, mult=mult))
        # opt = RSRFOptimizer(geom, trust_radius=0.5, update_trust=False,
                            # # hess_update="flowchart", max_micro_cycles=25,
                            # hess_update="flowchart", max_micro_cycles=1,
                            # max_cycles=150, thresh="gau")
        # opt.run()
        # converged = opt.is_converged
        # cycles = opt.cur_cycle + 1

        results[xyz_fn] = (converged, cycles)
        pprint(results)

        if converged:
            tot_cycles += cycles
        else:
            fails += 1
        print()
        if os.path.isfile("stop"):
            print("Found stop file")
            os.remove("stop")
            break
        sys.stdout.flush()
    print()

    pprint(results)
    print("Total cycles", tot_cycles)
    print("Fails", fails)


if __name__ == "__main__":
    run()
