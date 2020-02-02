#!/usr/bin/env python3

# [1] https://doi.org/10.1007/s002140050387
#     Bofill, 1998


import numpy as np

from pysisyphus.tsoptimizers.TSHessianOptimizer import TSHessianOptimizer


class RSIRFOptimizer(TSHessianOptimizer):

    def optimize(self):
        energy, gradient, H, eigvals, eigvecs = self.housekeeping()
        self.update_ts_mode(eigvals, eigvecs)

        self.log( "Using projection to construct image potential gradient "
                  f"and hessian for root {self.root}."
        )
        trans_vec = eigvecs[:,self.root]
        # Projection matrix to construct g* and H*
        P = np.eye(self.geometry.coords.size) - 2 * np.outer(trans_vec, trans_vec)
        H_star = P.dot(H)
        eigvals_, eigvecs_ = np.linalg.eigh(H_star)
        # Neglect small eigenvalues
        eigvals_, eigvecs_ = self.filter_small_eigvals(eigvals_, eigvecs_)

        # Transform gradient to basis of eigenvectors
        grad_star = eigvecs_.T.dot(P.dot(gradient))

        alpha = self.alpha0
        self.max_micro_cycles = 25
        for mu in range(self.max_micro_cycles):
            # assert alpha > 0, "alpha should not be negative"
            self.log(f"RS-IRFO micro cycle {mu:02d}, alpha={alpha:.6f}")
            H_aug = self.get_augmented_hessian(eigvals_, grad_star, alpha)
            rfo_step_, eigval_min, nu = self.solve_rfo(H_aug, "min")
            rfo_norm_ = np.linalg.norm(rfo_step_)
            self.log(f"norm(rfo step)={rfo_norm_:.6f}")

            if (rfo_norm_ < self.trust_radius) or abs(rfo_norm_ - self.trust_radius) <= 1e-3:
                step_ = rfo_step_
                break

            alpha_step = self.get_alpha_step(alpha, eigval_min, rfo_norm_, eigvals_, grad_star)
            alpha += alpha_step
            self.log("")

        # Transform back to original basis
        step = eigvecs_.dot(step_)

        quadratic_prediction = step @ gradient + 0.5 * step @ self.H @ step
        rfo_prediction = quadratic_prediction / (1 + step @ step)
        self.predicted_energy_changes.append(rfo_prediction)

        return step
