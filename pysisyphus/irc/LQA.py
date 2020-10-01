# [1] https://aip.scitation.org/doi/pdf/10.1063/1.459634?class=pdf
#     Page, 1990, Eq. 19 is missing a **2 after g'_0,i
# [2] https://aip.scitation.org/doi/10.1063/1.1724823
#     Hratchian, 2004

import numpy as np

from pysisyphus.optimizers.hessian_updates import bfgs_update
from pysisyphus.irc.IRC import IRC


class LQA(IRC):

    def __init__(self, geometry, N_euler=5000, **kwargs):
        super().__init__(geometry, **kwargs)

        self.N_euler = N_euler

    def step(self):
        mw_gradient = self.mw_gradient

        if len(self.irc_mw_gradients) > 1:
            dg = self.irc_mw_gradients[-1] - self.irc_mw_gradients[-2]
            dx = self.irc_mw_coords[-1] - self.irc_mw_coords[-2]
            dH, _ = bfgs_update(self.mw_hessian, dx, dg)
            self.mw_hessian += dH

        eigenvalues, eigenvectors = np.linalg.eigh(self.mw_hessian)
        # Drop small eigenvalues and corresponding eigenvectors
        small_vals = np.abs(eigenvalues) < 1e-8
        eigenvalues = eigenvalues[~small_vals]
        eigenvectors = eigenvectors[:,~small_vals]

        # t step for numerical integration
        dt = 1 / self.N_euler * self.step_length / np.linalg.norm(mw_gradient)

        # Transform gradient to eigensystem of the hessian
        mw_gradient_trans = eigenvectors.T @ mw_gradient

        t = dt
        cur_length = 0
        for i in range(self.N_euler):
            dsdt = np.sqrt(np.sum(mw_gradient_trans**2 * np.exp(-2*eigenvalues*t)))
            cur_length += dsdt * dt
            if cur_length > self.step_length:
                break
            t += dt
        alphas = (np.exp(-eigenvalues*t) - 1) / eigenvalues
        A = eigenvectors @ np.diag(alphas) @ eigenvectors.T
        step = A @ mw_gradient

        mw_coords = self.mw_coords.copy()
        self.mw_coords = mw_coords + step
