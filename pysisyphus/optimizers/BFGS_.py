#!/usr/bin/env python3

import numpy as np
# import matplotlib.pyplot as plt

from pysisyphus.helpers import procrustes, fit_rigid
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.line_search import line_search

np.set_printoptions(suppress=True, precision=4)

# [1] Nocedal, Wright - Numerical Optimization, 2006

class BFGS_(Optimizer):

    def __init__(self, geometry, alpha=1.0, **kwargs):
        super().__init__(geometry, alpha=alpha, **kwargs)

        self.eye = np.eye(len(self.geometry.coords))

    def prepare_opt(self):
        if self.is_cos and self.align:
            procrustes(self.geometry)
        self.H = self.eye

    def restrict_step(self, step):
        too_big = np.abs(step) > self.max_step
        self.log(f"Found {np.sum(too_big)} big step components.")
        signs = np.sign(step[too_big])
        step[too_big] = signs * self.max_step
        return step

    def optimize(self):
        if self.is_cos and self.align:
            _, rot_vec_lists, _ = fit_rigid(
                self.geometry,
                vector_lists=(self.coords, self.forces, self.steps)
            )
            self.coords, self.forces, self.steps = rot_vec_lists
        # realigned = False
        # if (self.cur_cycle > 0) and (self.cur_cycle % 10 == 0):
            # procrustes(self.geometry)
            # self.H = self.eye.copy()
            # print("realigned")
        forces = self.geometry.forces
        energy = self.geometry.energy

        self.forces.append(forces)
        self.energies.append(energy)

        # if self.cur_cycle > 0 and not realigned:
        if self.cur_cycle > 0:
            # coord_diff = self.coords[-1] - self.coords[-2]
            # grad_diff = self.forces[-2] - self.forces[-1]

            # curv_cond = coord_diff @  grad_diff
            # self.log(f"Curvature condition: {curv_cond:.6e}")

            # Scale initial hessian after first step and before first update
            # if self.cur_cycle == 1:
                # factor = (grad_diff @ coord_diff) / (grad_diff @ grad_diff)
                # self.H *= factor

            # BFGS update for the original hessian (NOT its inverse!)
            # first_term = ((self.H @ np.outer(coord_diff, coord_diff) @ self.H)
                          # / (coord_diff @ self.H @ coord_diff)
            # )
            # sec_term = np.outer(grad_diff, grad_diff) / (grad_diff @ coord_diff)
            # self.H = self.H - first_term + sec_term

            # BFGS update for the inverse hessian
            # rho = 1 / (grad_diff @ coord_diff)
            # A = self.eye - (rho * np.outer(coord_diff, grad_diff))
            # B = self.eye - (rho * np.outer(grad_diff, coord_diff))
            # C = rho * np.outer(coord_diff, coord_diff)
            # self.H = A @ self.H @ B + C
            coord_diffs = np.diff(self.coords, axis=0)
            grad_diffs = np.diff(-np.array(self.forces), axis=0)
            self.H = self.eye.copy()
            cdgd = list(zip(coord_diffs, grad_diffs))
            for coord_diff, grad_diff in cdgd[-10:]:
                rho = 1 / (grad_diff @ coord_diff)
                A = self.eye - (rho * np.outer(coord_diff, grad_diff))
                B = self.eye - (rho * np.outer(grad_diff, coord_diff))
                C = rho * np.outer(coord_diff, coord_diff)
                self.H = A @ self.H @ B + C
            
        # Step direction from the original Hessian
        # step_direction = np.linalg.pinv(self.H) @ forces

        # Step direction from the inverse Hessian
        step_direction = self.H @ forces
        alpha = self.alpha
        # alpha = line_search(f, df, x, p, f_0, df_0, alpha=1, alpha_max=1, c1=1e-4, c2=0.9)
        step = alpha * step_direction
        step = self.restrict_step(step)

        return step
