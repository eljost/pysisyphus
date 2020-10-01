#!/usr/bin/env python3

# [1] The Importance of Step Control in Optimization Methods

from math import sqrt

import numpy as np
from scipy.optimize import root_scalar

from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer


class RSAlgorithm(HessianOptimizer):

    def optimize(self):
        g = self.geometry.gradient
        self.forces.append(-self.geometry.gradient.copy())
        self.energies.append(self.geometry.energy)

        if self.cur_cycle > 0:
            self.update_trust_radius()
            self.update_hessian()

        H = self.H
        if self.geometry.internal:
            H = self.geometry.internal.project_hessian(H)

        vals, vecsT = np.linalg.eigh(H)
        # Exclude small eigenvalues and corresponding -vectors
        # small_inds = np.abs(vals) < 1e-6
        small_inds = np.abs(vals) < 1e-10
        self.log(f"Found {small_inds.sum()} small eigenvalues.")
        neg_inds = vals < 0
        neg = neg_inds.sum()

        # vals = vals[~small_inds]
        # vecsT = vecsT[:,~small_inds]
        g_ = vecsT.T.dot(g)
        self.log(f"shape(g_)={g.shape}, shape(vals)={vals.shape}, shape(vecsT)={vecsT.shape}")

        def get_step(lambda_):
            """Returns step for a given lambda"""
            # _ = -g_/(vals+lambda_)
            # TODO: remove offset
            _ = -g_/(vals+lambda_)
            return vecsT.dot(_)

        if neg == 0:
            self.log("Hessian is positive definite.")
            self.log("Checking pure QN-step ... ")
            step = get_step(0)
            step_norm = np.linalg.norm(step)
            if step_norm <= self.trust_radius:
                self.log(f"norm(step)={step_norm:.6f} is fine")
                predicted_change = step.dot(g) + 0.5 * step.dot(self.H).dot(step)
                self.predicted_energy_changes.append(predicted_change)
                return step
            self.log(f"norm(QN-step)={step_norm:.6f} is too big.")
        else:
            # Hessian is not positive definite
            self.log("Hessian is not positive definite.")
            smallest_eigval = vals[0]
            self.log(f"Smallest eigenvalue is {smallest_eigval:.6f}")

        def on_sphere_linear(lambda_):
            return 1/self.trust_radius - 1/np.linalg.norm(get_step(lambda_))

        _b1 = -vals[0]
        x0 = _b1 + 1e-3
        x1 = x0 + 1e-3

        # Defining a bracket using infinity (float("inf")) doesn't seem to work.
        # Instead we use a big number.
        upper_bracket = 1e10
        # Hessian  is positive definite
        if neg == 0:
            bracket = [0, upper_bracket]
        # Hessian has negative eigenvalues, is not positive definitie
        else:
            bracket = [_b1+1e-6, upper_bracket]
        sol = root_scalar(on_sphere_linear, x0=x0, x1=x1, bracket=bracket)
        if not sol.converged:
            raise Exception("Root search did not converge!")
        lambda_ = sol.root
        # Check shifted hessian for positive definiteness by adding the shift
        # to the smallest eigenvalue and see if this is > 0. If so we can use this
        # step (Case II).
        if vals[0] + lambda_ > 0:
            step = get_step(lambda_)
            step_norm = np.linalg.norm(step)
            self.log(f"Found valid step with λ={lambda_:.6f} and norm={step_norm:.6f}")
            predicted_change = step.dot(g) + 0.5 * step.dot(self.H).dot(step)
            self.predicted_energy_changes.append(predicted_change)
            return step

        import pdb; pdb.set_trace()
        self.log(f"Shifted hessian (λ={lambda_:.6f} is not positive definite!")
        self.log("Determining new step using second parameter τ.")
        # Shifted hessian is not positive definite (Case III).
        lambda_ = vals[0]
        # frac_sum = np.sum(g_[1:] / (vals[1:] - vals[0]))
        frac_sum = np.sum(g_[1:] / (vals[1:] - lambda_))
        tau = sqrt(self.trust_radius**2 - frac_sum**2)
        self.log(f"τ={tau:.6f}")

        # The second term is still in the eigensystem of the hessian, whereas the
        # first term is in the original space. So this won't work ...
        raise Exception("first term is in original space, second term is "\
                        "is still in eigenspace of hessian.")
        step = get_step(lambda_) + tau*g_[:,0]

        predicted_change = step.dot(g) + 0.5 * step.dot(self.H).dot(step)
        self.predicted_energy_changes.append(predicted_change)

        return step
