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

    def __init__(self, geometry, root=0, hessian_ref=None, prim_coord=None,
                 hessian_init="calc", hessian_update="bofill",
                 max_micro_cycles=50, trust_radius=0.3, **kwargs):

        # assert hessian_init == "calc", \
            # "TS-optimization should be started from a calculated hessian " \
            # "(hessian_init=\"calc\")!"
        assert hessian_update == "bofill", \
            "Bofill update is recommended in a TS-optimization."

        super().__init__(geometry, hessian_init=hessian_init,
                         hessian_update=hessian_update, **kwargs)

        self.root = int(root)
        self.hessian_ref = hessian_ref
        try:
            self.hessian_ref = np.loadtxt(self.hessian_ref)
            _ = self.geometry.coords.size
            expected_shape = (_, _)
            shape = self.hessian_ref.shape
            assert shape == expected_shape, \
                f"Shape of reference hessian {shape} doesn't match the expected "  \
                f"shape {expected_shape} of the hessian for the current coordinates!"
        except OSError as err:
            self.log(f"Tried to load reference hessian from '{self.hessian_ref}' "
                      "but the file could not be found.")
        except ValueError as err:
            self.log(f"No reference hessian provided.")
        self.prim_coord = prim_coord

        self.ts_mode = None
        self.max_micro_cycles = max_micro_cycles
        self.prim_contrib_thresh = 0.05

        self.alpha0 = 1

    def prepare_opt(self):
        super().prepare_opt()

        # Determiniation of initial mode either by using a provided
        # reference hessian, or by using a supplied root.

        eigvals, eigvecs = np.linalg.eigh(self.H)

        self.log_negative_eigenvalues(eigvals, "Initial ")

        self.log("Determining initial TS mode to follow uphill.")
        # Select an initial TS-mode by highest overlap with eigenvectors from
        # reference hessian.
        if self.prim_coord:
            prim_ind = self.geometry.internal.get_index_of_prim_coord(self.prim_coord)
            if prim_ind is None:
                raise Exception(f"Primitive internal {self.prim_coord} is not defined!")
            # Select row of eigenvector-matrix that corresponds to this coordinate
            prim_row = eigvecs[prim_ind]
            max_contrib_ind = np.abs(prim_row).argmax()
            self.root = max_contrib_ind

            big_contribs = np.abs(prim_row) > self.prim_contrib_thresh
            big_inds = np.arange(prim_row.size)[big_contribs]
            contrib_str = "\n".join(
                [f"\t{ind:03d}: {contrib:.4f}"
                 for ind, contrib in zip(big_inds, np.abs(prim_row)[big_contribs])]
            )

            self.log( "Highest absolute contribution of primitive internal coordinate "
                     f"{self.prim_coord} in mode {self.root:02d}.")
            self.log(f"Absolute contributions > {self.prim_contrib_thresh:.04f}"
                     f"\n{contrib_str}")
            used_str = f"contribution of primitive coordinate {self.prim_coord}"
        elif self.hessian_ref is not None:
            eigvals_ref, eigvecs_ref = np.linalg.eigh(self.hessian_ref)
            self.log_negative_eigenvalues(eigvals_ref, "Reference ")
            assert eigvals_ref[0] < -self.small_eigval_thresh
            ref_mode = eigvecs_ref[:,0]
            overlaps = np.einsum("ij,j->i", eigvecs.T, ref_mode)
            ovlp_str = np.array2string(overlaps, precision=4)
            self.log( "Overlaps between eigenvectors of current hessian "
                     f"TS mode from reference hessian:")
            self.log(f"\t{ovlp_str}")

            self.root = np.abs(overlaps).argmax()
            print( "Highest overlap between reference TS mode and "
                  f"eigenvector {self.root:02d}.")
            used_str = "overlap with reference TS mode"
        else:
            used_str = f"root={self.root}"
        self.log(f"Used {used_str} to select inital TS mode.")

        # Check if the selected mode (root) is a sensible choice.
        #
        # small_eigval_thresh is positive and we dont take the absolute value
        # of the eigenvalues. So we multiply small_eigval_thresh to get a
        # negative number.
        assert eigvals[self.root] < -self.small_eigval_thresh, \
             "Expected negative eigenvalue! Eigenvalue of selected TS-mode " \
            f"{self.root:02d} is above the the threshold of " \
            f"{self.small_eigval_thresh:.6e}!"

        # Select an initial TS-mode by root index. self.root may have been
        # modified by using a reference hessian.
        self.ts_mode = eigvecs[:,self.root]
        self.log(f"Using mode {self.root:02d} with eigenvalue "
                 f"{eigvals[self.root]:.6f} as TS mode.")
        self.log("")

    def update_ts_mode(self, eigvals, eigvecs):
        neg_eigval_inds = eigvals < -self.small_eigval_thresh
        neg_num = neg_eigval_inds.sum()
        assert neg_num >= 1, \
            "Need at least 1 negative eigenvalue for TS optimization."
        self.log_negative_eigenvalues(eigvals)

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
            # Symmetrize hessian, as the projection may break it?!
            H = (H_proj + H_proj.T) / 2

        eigvals, eigvecs = np.linalg.eigh(H)

        if self.geometry.coord_type == "cart":
            # Don't use eigenvectors that belong to very small eigenvalues,
            # as they probably belong to overall translations/rotations of
            # the molecule and may mess up the stepsize, by producing large steps.
            eigvals, eigvecs = self.filter_small_eigvals(eigvals, eigvecs)

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

        # Eq. (6) from [4] seems erronous ... the prediction is usually only ~50%
        # of the actual change ...
        # predicted_energy_change = 1/2 * (eigval_max / nu_max**2 + eigval_min / nu_min**2)
        # self.predicted_energy_changes.append(predicted_energy_change)

        quadratic_prediction = step @ gradient + 0.5 * step @ self.H @ step
        rfo_prediction = quadratic_prediction / (1 + step @ step)
        self.predicted_energy_changes.append(rfo_prediction)

        self.log("")
        return step
