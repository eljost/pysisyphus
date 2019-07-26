#!/usr/bin/env python3

# See [1] https://pubs.acs.org/doi/pdf/10.1021/j100247a015
#         Banerjee, 1985
#     [2] https://aip.scitation.org/doi/abs/10.1063/1.2104507
#         Heyden, 2005
#     [3] https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.540070402
#         Baker, 1985


import numpy as np

from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer


class PRFOptimizer(HessianOptimizer):
    """Optimizer to find first-order saddle points."""

    rfo_dict = {
        "min": (0, "min"),
        "max": (-1, "max"),
    }

    def __init__(self, geometry, root=0,
                 hessian_init="calc", hessian_update="bofill",
                 neg_eigval_thresh=-1e-8, **kwargs):

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
        eigenvalues, eigenvectors = np.linalg.eigh(rfo_mat)
        ind, verbose = self.rfo_dict[kind]
        # Eigenvalues and -values are sorted, so either use the first
        # (for minimization) or the last (for maximization) eigenvalue
        # and eigenvector.
        step = eigenvectors.T[ind]
        nu = step[-1]
        self.log(f"nu_{verbose}={nu:.4e}")
        # Scale eigenvectors corresponding to the largest (maximization)
        # or smallest (minimization) eigenvalue, so the last entry is 1.
        # The step then corresponds to the scaled eigenvector, without
        # the last element.
        step = step[:-1] / nu
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
            self.log(f"{i:02d}: {ovlp:.6f}")
        max_ovlp_ind = np.argmax(ovlps)
        max_ovlp = ovlps[max_ovlp_ind]
        self.log(f"Highest overlap: {max_ovlp:.6f}, mode {max_ovlp_ind}")
        self.log(f"Continuing with mode {max_ovlp_ind} as TS mode.")
        self.root = max_ovlp_ind
        self.ts_mode = eigvecs.T[max_ovlp_ind]

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

        step_max, eigval_max, nu_max = self.solve_rfo(max_mat, "max")
        step_max = step_max[0]
        step_min, eigval_min, nu_min = self.solve_rfo(min_mat, "min")

        # Assemble step from step_max and step_min
        step = np.zeros_like(gradient_trans)
        step[self.root] = step_max
        step[min_indices] = step_min
        # Right now the step is still given in the Hessians eigensystem. We
        # transform it back now.
        step = eigvecs.dot(step)
        step_norm = np.linalg.norm(step)

        # Restrict step_norm to the current trust radius
        if step_norm > self.trust_radius:
            step = step / step_norm * self.trust_radius
        self.log(f"norm(step)={np.linalg.norm(step):.6f}")

        predicted_change = step.dot(gradient) + 0.5 * step.dot(self.H).dot(step)
        self.predicted_energy_changes.append(predicted_change)

        self.log("")
        return step
