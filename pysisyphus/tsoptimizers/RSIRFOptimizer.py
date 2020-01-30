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

        # Augmented hessian
        dim_ = eigvals_.size + 1
        H_aug = np.zeros((dim_, dim_))
        H_aug[:dim_-1,:dim_-1] = np.diag(eigvals_)
        H_aug[-1,:-1] = grad_star
        H_aug[:-1,-1] = grad_star
        alpha = self.alpha0

        diag_indices = np.diag_indices(eigvals_.size)
        self.max_micro_cycles = 25
        for mu in range(self.max_micro_cycles):
            # assert alpha > 0, "alpha should not be negative"
            self.log(f"RS-IRFO micro cycle {mu:02d}, alpha={alpha:.6f}")
            # We only have to update one eigenvalue
            H_aug_scaled = H_aug.copy()
            H_aug_scaled[diag_indices] /= alpha
            H_aug_scaled[:-1,-1] /= alpha
            # import pdb; pdb.set_trace()
            rfo_step_, eigval_min, nu = self.solve_rfo(H_aug_scaled, "min")
            rfo_norm_ = np.linalg.norm(rfo_step_)
            self.log(f"norm(rfo step)={rfo_norm_:.6f}")

            if (rfo_norm_ < self.trust_radius) or abs(rfo_norm_ - self.trust_radius) <= 1e-3:
                step_ = rfo_step_
                break

            # Derivative of the squared step w.r.t. alpha
            numer = grad_star**2
            denom = (eigvals_ - eigval_min * alpha)**3
            quot = np.sum(numer / denom)
            self.log(f"quot={quot:.6f}")
            dstep2_dalpha = (2*eigval_min/(1+rfo_norm_**2 * alpha)
                             * np.sum(grad_star**2
                                      / ((eigvals_ - eigval_min * alpha)**3)
                               )
            )
            self.log(f"analytic deriv.={dstep2_dalpha:.6f}")
            # Update alpha
            alpha_step = (2*(self.trust_radius*rfo_norm_ - rfo_norm_**2)
                          / dstep2_dalpha
            )
            self.log(f"alpha_step={alpha_step:.4f}")
            alpha += alpha_step
            self.log("")

        # Transform back to original basis
        step = eigvecs_.dot(step_)

        quadratic_prediction = step @ gradient + 0.5 * step @ self.H @ step
        rfo_prediction = quadratic_prediction / (1 + step @ step)
        self.predicted_energy_changes.append(rfo_prediction)

        return step
