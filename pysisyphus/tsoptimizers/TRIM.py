#!/usr/bin/env python3

# [1] https://doi.org/10.1016/0009-2614(91)90115-P
#     Helgaker, 1991


import numpy as np

from scipy.optimize import newton
from pysisyphus.tsoptimizers.TSHessianOptimizer import TSHessianOptimizer


class TRIM(TSHessianOptimizer):

    def optimize(self):
        energy, gradient, H, eigvals, eigvecs = self.housekeeping()
        self.update_ts_mode(eigvals, eigvecs)

        self.log(f"Signs of eigenvalue and -vector of root {self.root} "
                  "will be reversed!")
        # Transform gradient to basis of eigenvectors
        gradient_ = eigvecs.T.dot(gradient)

        # Construct image function by inverting the signs of the eigenvalue and
        # -vector of the mode to follow uphill.
        eigvals_ = eigvals.copy()
        eigvals_[self.root] *= -1
        gradient_ = gradient_.copy()
        gradient_[self.root] *= -1

        def get_step(mu):
            zetas = -gradient_ / (eigvals_ - mu)
            # Replace nan with 0.
            zetas = np.nan_to_num(zetas)
            # Transform to original basis
            step = eigvecs * zetas
            step = step.sum(axis=1)
            return step

        def get_step_norm(mu):
            return np.linalg.norm(get_step(mu))

        def func(mu):
            return get_step_norm(mu) - self.trust_radius

        mu = 0
        norm0 = get_step_norm(mu)
        if norm0 > self.trust_radius:
            mu, res = newton(func, x0=mu, full_output=True)
            assert res.converged
            self.log(f"Using levelshift of μ={mu:.4f}")
        else:
            self.log("Took pure newton step without levelshift")

        step = get_step(mu)
        step_norm = np.linalg.norm(step)
        self.log(f"norm(step)={step_norm:.6f}")

        quadratic_prediction = step @ gradient + 0.5 * step @ self.H @ step
        self.predicted_energy_changes.append(quadratic_prediction)

        return step
