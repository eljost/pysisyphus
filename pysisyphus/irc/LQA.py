#!/usr/bin/env python3

# [1] https://aip.scitation.org/doi/pdf/10.1063/1.459634?class=pdf
#     Page, 1990, Eq. 19 is missing a **2 after g'_0,i
# [2] https://aip.scitation.org/doi/10.1063/1.1724823
#     Hratchian, 2004

import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.hessian_updates import bfgs_update
from pysisyphus.irc.IRC import IRC


class LQA(IRC):

    def __init__(self, geometry, N_euler=5000, **kwargs):
        super().__init__(geometry, **kwargs)

        self.N_euler = N_euler

    def prepare(self, direction):
        super().prepare(direction)

        mm_sqrt_inv = self.geometry.mm_sqrt_inv
        self.mw_hessian = mm_sqrt_inv @ self.hessian @ mm_sqrt_inv
        self.mw_gradients = list()

    def step(self):
        gradient = self.mw_gradient
        self.mw_gradients.append(gradient)
        # hessian = self.geometry.mw_hessian
        hessian = self.mw_hessian

        if len(self.mw_gradients) > 1:
            dg = self.mw_gradients[-1] - self.mw_gradients[-2]
            dx = self.irc_mw_coords[-1] - self.irc_mw_coords[-2]
            dH, _ = bfgs_update(self.mw_hessian, dx, dg)
            self.mw_hessian += dH

        eigenvalues, eigenvectors = np.linalg.eigh(hessian)
        # Drop small eigenvalues and corresponding eigenvectors
        small_vals = np.abs(eigenvalues) < 1e-8
        eigenvalues = eigenvalues[~small_vals]
        eigenvectors = eigenvectors[:,~small_vals]

        # t step for numerical integration
        dt = 1 / self.N_euler * self.step_length / np.linalg.norm(gradient)

        # Transform gradient to eigensystem of the hessian
        gradient_trans = eigenvectors.T @ gradient

        t = dt
        _ = 0
        for i in range(self.N_euler):
            dsdt = np.sqrt(np.sum(gradient_trans**2 * np.exp(-2*eigenvalues*t)))
            _ += dsdt * dt
            if _ > self.step_length:
                break
            t += dt
        alphas = (np.exp(-eigenvalues*t) - 1) / eigenvalues
        A = eigenvectors @ np.diag(alphas) @ eigenvectors.T
        step = A @ gradient
        norm = np.linalg.norm(step)

        mw_coords = self.mw_coords.copy()
        self.mw_coords = mw_coords + step
