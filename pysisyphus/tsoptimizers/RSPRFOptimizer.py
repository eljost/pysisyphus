#!/usr/bin/env python3

# See [1] https://pubs.acs.org/doi/pdf/10.1021/j100247a015
#         Banerjee, 1985
#     [2] https://aip.scitation.org/doi/abs/10.1063/1.2104507
#         Heyden, 2005
#     [3] https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.540070402
#         Baker, 1985
#     [4] https://link.springer.com/article/10.1007/s002140050387
#         Besalu, 1998


import numpy as np

from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer


class RSPRFOptimizer(HessianOptimizer):
    """Optimizer to find first-order saddle points."""

    rfo_dict = {
        "min": (0, "min"),
        "max": (-1, "max"),
    }

    def __init__(self, geometry, root=0,
                 hessian_init="calc", hessian_update="bofill",
                 neg_eigval_thresh=-1e-8, max_micro_cycles=50, **kwargs):

        assert hessian_init == "calc", \
            "TS-optimization should be started from a calculated hessian " \
            "(hessian_init=\"calc\")!"
        assert hessian_update == "bofill", \
            "Bofill update is recommended in a TS-optimization."

        super().__init__(geometry, hessian_init=hessian_init,
                         hessian_update=hessian_update, **kwargs)

        self.root = int(root)
        self.ts_mode = None
        self.neg_eigval_thresh = float(neg_eigval_thresh)
        self.max_micro_cycles = max_micro_cycles

        self.alpha0 = 1

    def prepare_opt(self):
        super().prepare_opt()
        eigvals, eigvecs = np.linalg.eigh(self.H)
        # Check if the selected mode is a sensible choice
        assert eigvals[self.root] < self.neg_eigval_thresh, \
             "Expected negative eigenvalue! Eigenvalue of selected TS-mode " \
            f"{self.root:02d} is above the the threshold of " \
            f"{self.neg_eigval_thresh:.6e}!"
        # Select an initial TS-mode
        self.ts_mode = eigvecs[:,self.root]

    def solve_rfo(self, rfo_mat, kind="min"):
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

    def update_ts_mode(self, eigvals, eigvecs):
        neg_eigval_inds = eigvals < -1e-8
        neg_num = neg_eigval_inds.sum()
        assert neg_num >= 1, \
            "Need at least 1 negative eigenvalue for TS optimization."
        eigval_str = np.array2string(eigvals[neg_eigval_inds], precision=6)
        self.log(f"Found {neg_num} negative eigenvalue(s): {eigval_str}")
        # Select TS mode with biggest overlap to the previous TS mode
        self.log("Overlaps of previous TS mode with current imaginary mode(s):")
        ovlps = [np.abs(imag_mode.dot(self.ts_mode)) for imag_mode in eigvecs.T[:neg_num]]
        for i, ovlp in enumerate(ovlps):
            self.log(f"\t{i:02d}: {ovlp:.6f}")
        max_ovlp_ind = np.argmax(ovlps)
        max_ovlp = ovlps[max_ovlp_ind]
        self.log(f"Highest overlap: {max_ovlp:.6f}, mode {max_ovlp_ind}")
        self.log(f"Continuing with mode {max_ovlp_ind} as TS mode.")
        self.root = max_ovlp_ind
        self.ts_mode = eigvecs.T[max_ovlp_ind]

    # def optimize(self):
        # gradient = self.geometry.gradient
        # self.forces.append(-self.geometry.gradient)
        # self.energies.append(self.geometry.energy)

        # if self.cur_cycle > 0:
            # self.update_trust_radius()
            # self.update_hessian()

        # H = self.H
        # if self.geometry.coord_type == "redund":
            # H_proj = self.geometry.internal.project_hessian(self.H)
            # # Symmetrize hessian, as the projection probably breaks it?!
            # H = (H_proj + H_proj.T) / 2

        # eigvals, eigvecs = np.linalg.eigh(H)
        # self.update_ts_mode(eigvals, eigvecs)

        # # Transform gradient to eigensystem of hessian
        # gradient_trans = eigvecs.T.dot(gradient)

        # # Minimize energy along all modes, except the TS-mode
        # min_indices = [i for i in range(gradient_trans.size) if i != self.root]
        # min_mat = np.asarray(np.bmat((
            # (np.diag(eigvals[min_indices]), gradient_trans[min_indices,None]),
            # (gradient_trans[None,min_indices], [[0]])
        # )))
        # # Maximize energy along the chosen TS mode. The matrix is hardcoded
        # # as 2x2, so only first-order saddle point searches are supported.
        # max_mat = np.array(((eigvals[self.root], gradient_trans[self.root]),
                           # (gradient_trans[self.root], 0)))

        # step_max, eigval_max, nu_max = self.solve_rfo(max_mat, "max")
        # step_max = step_max[0]
        # step_min, eigval_min, nu_min = self.solve_rfo(min_mat, "min")

        # # Assemble step from step_max and step_min
        # step = np.zeros_like(gradient_trans)
        # step[self.root] = step_max
        # step[min_indices] = step_min
        # # Right now the step is still given in the Hessians eigensystem. We
        # # transform it back now.
        # step = eigvecs.dot(step)
        # step_norm = np.linalg.norm(step)

        # # Restrict step_norm to the current trust radius
        # if step_norm > self.trust_radius:
            # step = step / step_norm * self.trust_radius
        # self.log(f"norm(step)={np.linalg.norm(step):.6f}")

        # predicted_change = step.dot(gradient) + 0.5 * step.dot(self.H).dot(step)
        # self.predicted_energy_changes.append(predicted_change)

        # self.log("")
        # return step

    def optimize(self):
        gradient = self.geometry.gradient
        self.forces.append(-self.geometry.gradient)
        self.energies.append(self.geometry.energy)

        if self.cur_cycle > 0:
            self.update_trust_radius()
            self.update_hessian()

        H = self.H
        if self.geometry.internal:
            H_proj = self.geometry.internal.project_hessian(self.H)
            # Symmetrize hessian, as the projection probably breaks it?!
            H = (H_proj + H_proj.T) / 2

        eigvals, eigvecs = np.linalg.eigh(H)

        if self.geometry.coord_type == "cart":
            # Poor mans Eckart projection ... or how to neglect translation
            # and rotation.
            #
            # Don't use eigenvectors that belong to very small eigenvalues,
            # as they belong to overall translations/rotations of the molecule
            small_inds = np.abs(eigvals) < 1e-8
            eigvals = eigvals[~small_inds]
            eigvecs = eigvecs[:,~small_inds]
            small_num = sum(small_inds)
            assert small_num in (5, 6), \
                 "Expected 5 or 6 small eigenvalues in cartesian hessian " \
                f"but found {small_num}!"

        # Calculate overlaps between (updated) hessian and current TS-mode to
        # determine new TS-mode.
        self.update_ts_mode(eigvals, eigvecs)

        # Transform gradient to eigensystem of hessian
        gradient_trans = eigvecs.T.dot(gradient)

        # Minimize energy along all modes, except the TS-mode
        min_indices = [i for i in range(gradient_trans.size) if i != self.root]
        min_mat = np.asarray(np.bmat((
            (np.diag(eigvals[min_indices]), gradient_trans[min_indices,None]),
            (gradient_trans[None,min_indices], [[0]])
        )))
        # Maximize energy along the chosen TS mode. The matrix is hardcoded
        # as 2x2, so only first-order saddle point searches are supported.
        max_mat = np.array(((eigvals[self.root], gradient_trans[self.root]),
                           (gradient_trans[self.root], 0)))

        min_diag_indices = np.diag_indices(min_mat.shape[0])

        """In the RS-(P)RFO method we have to scale the matrices with alpha.
        Unscaled matrix (Eq. 8) in [1]:
            (H  g) (x)          (S 0) (x)
                       = lambda
            (g+ 0) (1)          (0 1) (1)
        with
            S = alpha * Identity matrix
        and multiplying from the left with the inverse of the scaling matrix
            (1/alpha 0)

            (0       1)
        we get
            (1/alpha 0) (H  g) (x)          (x)
                                   = lambda
            (0       1) (g+ 0) (1)          (1)
        eventually leading to the scaled matrix:
            (H/alpha  g/alpha) (x)          (x)
                                   = lambda     .
            (g+             0) (1)          (1)
        """
        alpha = self.alpha0
        for mu in range(self.max_micro_cycles):
            self.log(f"RS-PRFO micro cycle {mu:02d}, alpha={alpha:.6f}")
            max_mat_scaled = max_mat.copy()
            # Create the scaled max. matrix
            max_mat_scaled[0, 0] /= alpha
            max_mat_scaled[0, 1] /= alpha
            step_max, eigval_max, nu_max = self.solve_rfo(max_mat_scaled, "max")
            step_max = step_max[0]
            # Create the scaled min. matrix
            scaled_min_eigenvalues = np.zeros_like(eigvals)
            scaled_min_eigenvalues[:-1] = eigvals[min_indices] / alpha
            min_mat_scaled = min_mat.copy()
            min_mat_scaled[min_diag_indices] = scaled_min_eigenvalues
            min_mat_scaled[:-1,-1] /= alpha
            step_min, eigval_min, nu_min = self.solve_rfo(min_mat_scaled, "min")

            # As of Eq. (8a) of [4] max_eigval and min_eigval also
            # correspond to:
            # max_eigval = -forces_trans[self.root] * max_step
            # min_eigval = -forces_trans[min_indices].dot(min_step)

            # Create the full PRFO step
            step = np.zeros_like(gradient_trans)
            step[self.root] = step_max
            step[min_indices] = step_min
            step_norm = np.linalg.norm(step)

            inside_trust = step_norm <= self.trust_radius
            if inside_trust:
                self.log("Restricted step satisfied the trust radius.")
                self.log(f"Micro-cycles converged in cycle {mu:02d} with "
                         f"alpha={alpha:.6f}!")
                break

            # Derivative of the squared step w.r.t. alpha
            # max subspace
            dstep2_dalpha_max = (2*eigval_max/(1+step_max**2 * alpha)
                                 * gradient_trans[self.root]**2
                                 / (eigvals[self.root] - eigval_max * alpha)**3
            )
            # min subspace
            dstep2_dalpha_min = (2*eigval_min/(1+step_min.dot(step_min) * alpha)
                                 * np.sum(gradient_trans[min_indices]**2
                                          / (eigvals[min_indices] - eigval_min * alpha)**3
                                 )
            )
            dstep2_dalpha = dstep2_dalpha_max + dstep2_dalpha_min
            # Update alpha
            alpha_step = (2*(self.trust_radius*step_norm - step_norm**2)
                          / dstep2_dalpha
            )
            alpha += alpha_step
            assert alpha > 0, "alpha must not be negative!"

        # Right now the step is still given in the Hessians eigensystem. We
        # transform it back now.
        step = eigvecs.dot(step)
        step_norm = np.linalg.norm(step)

        # With max_micro_cycles = 1 the RS part is disabled and the step
        # probably isn't scaled correctly in the one micro cycle.
        # In this case we use a naive scaling if the step is too big.
        if (self.max_micro_cycles == 1) and (step_norm > self.trust_radius):
            step = step / step_norm * self.trust_radius
        self.log(f"norm(step)={np.linalg.norm(step):.6f}")

        predicted_energy_change = 1/2 * eigval_max / nu_max**2 + eigval_min / nu_min**2
        self.predicted_energy_changes.append(predicted_energy_change)

        self.log("")
        return step
