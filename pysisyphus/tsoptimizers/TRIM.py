#!/usr/bin/env python3

import numpy as np

from scipy.optimize import newton
from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer


class TRIM(HessianOptimizer):

    def optimize(self):
        H = self.geometry.hessian
        gradient = self.geometry.gradient
        self.forces.append(-gradient)
        self.energies.append(self.geometry.energy)
        eigvals, eigvecs = np.linalg.eigh(H)

        root = 0
        assert root == 0
        self.log(f"Signs of eigenvalue and -vector of root {root} "
                  "will be reversed!")

        # Transform gradient to basis of eigenvectors
        gradient_ = eigvecs.T.dot(gradient)

        # Construct image function by inverting the signs of the eigenvalue and
        # -vector of the mode to follow uphill.
        eigvals_ = eigvals.copy()
        eigvals_[root] *= -1
        gradient_ = gradient_.copy()
        gradient_[root] *= -1

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
        return step
