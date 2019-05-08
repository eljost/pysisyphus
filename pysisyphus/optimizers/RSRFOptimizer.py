#!/usr/bin/env python3

# See [1] https://pubs.acs.org/doi/pdf/10.1021/j100247a015
#         Banerjee, 1985
#     [2] https://aip.scitation.org/doi/abs/10.1063/1.2104507
#         Heyden, 2005
#     [3] https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.540070402
#         Baker, 1985
#     [4] 10.1007/s002140050387
#         Bofill, 1998, Restricted-Step-RFO
#     [5] https://link.springer.com/article/10.1007/s00214-016-1847-3
#         Birkholz, 2016

import numpy as np

from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer


class RSRFOptimizer(HessianOptimizer):
    """Optimizer to find first-order saddle points."""

    rfo_dict = {
        "min": (0, "min"),
        "max": (-1, "max"),
    }

    def __init__(self, geometry, max_micro_cycles=50, **kwargs):
        super().__init__(geometry, **kwargs)

        self.max_micro_cycles = int(max_micro_cycles)
        assert max_micro_cycles >= 1

        self.alpha0 = 1
        self.alpha_max = 1e8

    def solve_rfo(self, rfo_mat, kind="min"):
        # So if I use eig instead of eigh here it even works ...
        # my bad, ahhh! The unscaled RFO matrix may be symmetric,
        # but the scaled ones aren't anymore.
        eigenvalues, eigenvectors = np.linalg.eig(rfo_mat)
        eigenvalues = eigenvalues.real
        eigenvectors = eigenvectors.real
        sorted_inds = np.argsort(eigenvalues)

        # Depending on wether we want to minimize (maximize) along
        # the mode(s) in the rfo mat we have to select the smallest
        # (biggest) eigenvalue and corresponding eigenvector.
        first_or_last, verbose = self.rfo_dict[kind]
        ind = sorted_inds[first_or_last]
        # Given sorted eigenvalue-indices (sorted_inds) use the first
        # (smallest eigenvalue) or the last (largest eigenvalue) index.
        step_nu = eigenvectors.T[ind]
        nu = step_nu[-1]
        self.log(f"nu_{verbose}={nu:.4e}")
        # Scale eigenvector so that its last element equals 1. The
        # final is step is the scaled eigenvector without the last element.
        step = step_nu[:-1] / nu
        eigval = eigenvalues[ind]
        self.log(f"eigenvalue_{verbose}={eigval:.4e}")
        return step, eigval, nu

    def optimize(self):
        forces = self.geometry.forces
        self.forces.append(forces)
        self.energies.append(self.geometry.energy)

        if self.cur_cycle > 0:
            self.update_trust_radius()
            self.update_hessian()

        H = self.H
        if self.geometry.internal:
            H = self.geometry.internal.project_hessian(self.H)
        eigvals, eigvecs = np.linalg.eigh(H)

        # Transform to eigensystem of hessian
        forces_trans = eigvecs.T.dot(forces)

        # Minimize energy along all modes
        min_mat = np.asarray(np.bmat((
            (np.diag(eigvals), -forces_trans[:,None]),
            (-forces_trans[None,:], [[0]])
        )))

        alpha = self.alpha0
        min_diag_indices = np.diag_indices(eigvals.size)
        for mu in range(self.max_micro_cycles):
            assert alpha > 0, "alpha should not be negative"
            self.log(f"RS-RFO micro cycle {mu:02d}, alpha={alpha:.6f}")
            # We only have to update one eigenvalue
            min_mat_scaled = min_mat.copy()
            min_mat_scaled[min_diag_indices] /= alpha
            min_mat_scaled[:-1,-1] /= alpha
            rfo_step, eigval_min, nu_min = self.solve_rfo(min_mat_scaled, "min")

            # As of Eq. (8a) of [4] max_eigval and min_eigval also
            # correspond to:
            # eigval_min_ = -forces_trans.dot(rfo_step)
            # np.testing.assert_allclose(eigval_min, eigval_min_)

            # Create the full PRFO step
            rfo_norm = np.linalg.norm(rfo_step)
            self.log(f"rfo_norm={rfo_norm:.6f}")

            inside_trust = rfo_norm < self.trust_radius + 1e-3
            if inside_trust:
                self.log("step is inside trust radius. breaking.")
                break
            elif alpha > self.alpha_max:
                print("alpha > alpha_max. breaking.")
                break

            # Derivative of the squared step w.r.t. alpha
            tval = 2*eigval_min/(1+rfo_norm**2 * alpha)
            numer = forces_trans**2
            denom = (eigvals - eigval_min * alpha)**3
            quot = np.sum(numer / denom)
            self.log(f"quot={quot:.6f}")
            dstep2_dalpha = (2*eigval_min/(1+rfo_norm**2 * alpha)
                             * np.sum(forces_trans**2
                                      / ((eigvals - eigval_min * alpha)**3)
                               )
            )
            self.log(f"analytic deriv.={dstep2_dalpha:.6f}")
            # Update alpha
            alpha_step = (2*(self.trust_radius*rfo_norm - rfo_norm**2)
                          / dstep2_dalpha
            )
            self.log(f"alpha_step={alpha_step:.4f}")
            alpha += alpha_step
            self.log("")

        # Right now the step is still given in the Hessians eigensystem. We
        # transform it back now.
        step = eigvecs.dot(rfo_step)
        step_norm = np.linalg.norm(step)
        # This would correspond to "pure" RFO without the iterative
        # step-restriction. Here we will just scale down the step, if it
        # is too big.
        if self.max_micro_cycles == 1 and step_norm > self.trust_radius:
            self.log("Scaled down step")
            step = step / step_norm * self.trust_radius
            step_norm = np.linalg.norm(step)

        self.log(f"norm(step)={np.linalg.norm(step):.6f}")

        # Calculating the energy change from eigval_min and nu_min seems to give
        # big problems.
        # predicted_energy_change = 1/2 * eigval_min / nu_min**2
        predicted_change = step.dot(-forces) + 0.5 * step.dot(self.H).dot(step)
        self.predicted_energy_changes.append(predicted_change)

        self.log("")
        return step
