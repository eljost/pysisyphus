#!/usr/bin/env python3

import numpy as np

from pysisyphus.helpers import geom_from_library
from pysisyphus.calculators.XTB import XTB
from pysisyphus.optimizers.hessian_updates import bfgs_update


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


def run():
    fn = "codein.xyz"
    geom = geom_from_library(fn, coord_type="redund")
    geom.set_calculator(XTB(pal=4))
    max_cycles = 50
    grads = list()
    steps = list()
    H = geom.get_initial_hessian()
    for i in range(max_cycles):
        grad = -geom.forces
        grads.append(grad)

        norm = np.linalg.norm(grad)
        if norm < 1e-4:
            print("Converged")
            break

        if i > 0:
            dx = steps[-1]
            dg = grads[-1] - grads[-2]
            dH, _ = bfgs_update(H, dx, dg)
            H += dH
        H_proj = H.copy()
        H_proj = geom.internal.project_hessian(H_proj)
        step = rfo(grad, H_proj, trust=0.3)
        steps.append(step)
        max_g = np.abs(grad).max()
        rms_g = np.sqrt(np.mean(grad**2))
        max_s = np.abs(step).max()
        rms_s = np.sqrt(np.mean(step**2))
        print(f"{i:02d}: max(f)={max_g:.6f}, rms(f)={rms_g:.6f}, "
              f"max(step)={max_s:.6f}, rms(step)={rms_s:.6f}")
        new_coords = geom.coords + step
        geom.coords = new_coords


if __name__ == "__main__":
    run()
