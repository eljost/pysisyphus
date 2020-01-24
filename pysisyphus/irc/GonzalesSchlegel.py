#!/usr/bin/env python3

import copy

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import newton

from pysisyphus.Geometry import Geometry
from pysisyphus.irc.IRC import IRC
from pysisyphus.optimizers.hessian_updates import bfgs_update
from pysisyphus.TableFormatter import TableFormatter

# [1] An improved algorithm for reaction path following
# http://aip.scitation.org/doi/pdf/10.1063/1.456010
# [2] Extension to internal coordinates (not implemented)
# https://pubs.acs.org/doi/pdf/10.1021/j100377a021

class GonzalesSchlegel(IRC):

    def __init__(self, geometry, max_micro_steps=20, **kwargs):
        super().__init__(geometry, **kwargs)

        self.max_micro_steps = max_micro_steps

        self.pivot_coords = list()
        self.micro_coords = list()

        micro_header = "# |dx| |tangent|".split()
        micro_fmts = ["d", ".2E", ".3E"]
        self.micro_formatter = TableFormatter(micro_header, micro_fmts, 10)

    def micro_step(self):
        """Constrained optimization on a hypersphere."""
        eye = np.eye(self.displacement.size)

        gradient = self.mw_gradient
        gradient_diff = gradient - self.prev_grad
        coords_diff = self.mw_coords - self.prev_coords
        self.prev_grad = gradient
        # Without copy we would only store the reference...
        self.prev_coords = self.mw_coords.copy()

        dH, _ = bfgs_update(self.hessian, coords_diff, gradient_diff)
        self.hessian += dH
        eigvals, eigvecs = np.linalg.eig(self.hessian)
        hessian_inv = np.linalg.pinv(self.hessian)

        def lambda_func(lambda_):
            # Eq. (11) in [1]
            # (H - λI)^-1
            hmlinv = np.linalg.pinv(self.hessian - eye*lambda_)
            # (g - λp)
            glp = gradient - self.displacement*lambda_
            tmp = self.displacement - hmlinv.dot(glp)
            return tmp.dot(tmp) - 0.25*(self.step_length**2)

        smallest_eigval = np.sort(eigvals)[0]
        # Initial guess for λ.
        # λ must be smaller then the smallest eigenvector
        lambda_ = np.sort(eigvals)[0]
        lambda_ *= 1.5 if (lambda_ < 0) else 0.5
        # Find the root with scipy
        lambda_ = newton(lambda_func, lambda_, maxiter=500)

        # Calculate dx from optimized lambda
        dx = -np.dot(
                np.linalg.inv(self.hessian-lambda_*eye),
                gradient-lambda_*self.displacement
        )
        self.displacement += dx
        self.mw_coords += dx

        displ_norm = np.linalg.norm(self.displacement)
        tangent = gradient - gradient.dot(self.displacement)/displ_norm * gradient

        return dx, tangent

    def step(self):
        grad0 = self.mw_gradient
        grad0_norm = np.linalg.norm(grad0)
        # For the BFGS update in the first micro step we use the original
        # point and the initial guess to calculate gradient and
        # coordinate differences.
        self.prev_grad = grad0
        self.prev_coords = self.mw_coords

        # Take a step against the gradient to the pivot point x*_k+1.
        pivot_step = -0.5*self.step_length * grad0/grad0_norm
        pivot_coords = self.mw_coords + pivot_step
        self.pivot_coords.append(pivot_coords)

        # Make initial guess for x'_k+1. Here we take another half
        # step from the pivot point.
        self.mw_coords = pivot_coords + pivot_step
        # Initial displacement p' from the pivot point
        self.displacement = pivot_step

        micro_coords_ = list()
        i = 0
        # self.table.print(f"Microiterations for step {self.cur_cycle}")
        self.table.print(self.micro_formatter.header)
        while True:
            if i == self.max_micro_steps:
                self.logger.warning("Max micro cycles exceeded!")
                break
            try:
                dx, tangent = self.micro_step()
            except RuntimeError:
                print("Constrained search did not converge!")
                self.converged = True
                return
            micro_coords_.append(self.mw_coords)
            norm_dx = np.linalg.norm(dx)
            norm_tangent = np.linalg.norm(tangent)
            self.table.print(self.micro_formatter.line(i+1, norm_dx, norm_tangent))

            if (np.linalg.norm(dx) <= 1e-3):
                break
            i += 1

        self.micro_coords.append(np.array(micro_coords_))

    def postprocess(self):
        self.pivot_coords = np.array(self.pivot_coords)
        self.micro_coords = np.array(self.micro_coords)
