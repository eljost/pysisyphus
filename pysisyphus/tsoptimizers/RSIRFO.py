#!/usr/bin/env python3

# [1] https://doi.org/10.1007/s002140050387
#     Bofill, 1998


import numpy as np

from scipy.optimize import newton
from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer


class RSIRFO(HessianOptimizer):
    # TODO: reuse methods from RSPRFO to select inital root etc.

    def prepare_opt(self):
        super().prepare_opt()
        self.I = np.eye(self.geometry.coords.size)

    def optimize(self):
        H = self.geometry.hessian
        gradient = self.geometry.gradient
        self.forces.append(-gradient)
        self.energies.append(self.geometry.energy)
        eigvals, eigvecs = np.linalg.eigh(H)
        # # Neglect small eigenvalues
        # eigvals, eigvecs = self.filter_small_eigvals(eigvals, eigvecs)

        root = 0
        assert root == 0
        trans_vec = eigvecs[:,root]
        # Projection matrix to construct g* and H*
        P = self.I - 2 * np.outer(trans_vec, trans_vec)
        H_star = P.dot(H)
        eigvals_, eigvecs_ = np.linalg.eigh(H_star)
        # Neglect small eigenvalues
        eigvals_, eigvecs_ = self.filter_small_eigvals(eigvals_, eigvecs_)

        # Transform gradient to basis of eigenvectors
        grad_star = eigvecs_.T.dot(P.dot(gradient))

        # Augmented hessian
        H_aug = np.array(np.bmat(
                            ((np.diag(eigvals_), grad_star[:, None]),
                             (grad_star[None, :], [[0]]))
        ))
        # alpha = self.alpha0
        alpha = 1

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
            tval = 2*eigval_min / (1+rfo_norm_**2 * alpha)
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

            # import pdb; pdb.set_trace()
            pass

        # Transform back to original basis
        step = eigvecs_.dot(step_)

        return step

        step_norm = np.linalg.norm(step)
        if step_norm > self.trust_radius:
            step = step / step_norm * self.trust_radius
        return step
