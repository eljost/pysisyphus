#!/usr/bin/env python3

# See [1] https://pubs.acs.org/doi/pdf/10.1021/j100247a015
#         Banerjee, 1985
#     [2] https://aip.scitation.org/doi/abs/10.1063/1.2104507
#         Heyden, 2005
#     [3] https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.540070402
#         Baker, 1985
#     [4] 10.1007/s002140050387
#         Bofill, 1998, Restricted-Step-RFO

import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer


class RSPRFOptimizer(Optimizer):
    """Optimizer to find first-order saddle points."""

    rfo_dict = {
        "min": (0, "min"),
        "max": (-1, "max"),
    }

    def __init__(self, geometry, root=0, trust_radius=.1, recalc_hess=None,
                 scaling_factor=2, max_micro_cycles=50,
                 **kwargs):
        super().__init__(geometry, **kwargs)

        self.root = int(root)
        self.trust_radius0 = trust_radius
        self.recalc_hess = recalc_hess
        self.scaling_factor = scaling_factor
        self.max_micro_cycles = max_micro_cycles

        self.alpha0 = 1
        self.trust_radius = self.trust_radius0
        self.H = None
        self.ts_mode = None
        self.predicted_energy_changes = list()

    def prepare_opt(self):
        self.H = self.geometry.hessian
        eigvals, eigvecs = np.linalg.eigh(self.H)
        self.ts_mode = eigvecs[:,self.root]

    def bofill_update(self, H, dx, dg):
        dgHdx = dg - H.dot(dx)

        # Symmetric, rank-one (SR1) update
        sr1 = np.outer(dgHdx, dgHdx) / dgHdx.dot(dx)

        # Powell update
        powell_1 = (np.outer(dgHdx, dx) + np.outer(dx, dgHdx)) / dx.dot(dx)
        powell_2 = dgHdx.dot(dx) * np.outer(dx, dx) / dx.dot(dx)**2
        powell = powell_1 - powell_2

        # Bofill mixing-factor
        mix = dgHdx.dot(dx)**2 / (dgHdx.dot(dgHdx) * dx.dot(dx))

        # Bofill update
        bofill_update = (mix * sr1) + (1 - mix)*(powell)

        return bofill_update

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

    def optimize(self):
        forces = self.geometry.forces
        self.forces.append(forces)
        self.energies.append(self.geometry.energy)

        if (self.recalc_hess and (self.cur_cycle > 1)
            and (self.cur_cycle % self.recalc_hess) == 0):
            self.log("Recalculating exact hessian")
            self.H = self.geometry.hessian
        elif len(self.coords) > 1:
            # Gradient difference
            dg = -(self.forces[-1] - self.forces[-2])
            # Coordinate difference
            dx = self.coords[-1] - self.coords[-2]
            self.H += self.bofill_update(self.H, dx, dg)
            self.log("Did Bofill hessian update.")
        if self.cur_cycle > 1:
            actual_energy_change = self.energies[-1] - self.energies[-2]
            predicted_energy_change = self.predicted_energy_changes[-1]
            coeff = actual_energy_change / predicted_energy_change
            print("coeff", coeff)

        eigvals, eigvecs = np.linalg.eigh(self.H)
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

        # Transform to eigensystem of hessian
        forces_trans = eigvecs.T.dot(forces)

        # Maximize energy along the chosen TS mode. The matrix is hardcoded
        # as 2x2, so only first-order saddle point searches are supported.
        max_mat = np.array(((eigvals[self.root], -forces_trans[self.root]),
                           (-forces_trans[self.root], 0)))
        # Minimize energy along all modes, except the TS mode
        min_indices = [i for i in range(forces.size) if i != self.root]
        min_mat = np.asarray(np.bmat((
            (np.diag(eigvals[min_indices]), -forces_trans[min_indices,None]),
            (-forces_trans[None,min_indices], [[0]])
        )))

        alpha = self.alpha0
        min_diag_indices = np.diag_indices(min_mat.shape[0])
        for mu in range(self.max_micro_cycles):
            self.log(f"RS-RFO micro cycle {mu:02d}")
            max_mat_scaled = max_mat.copy()
            # We only have to update one eigenvalue
            max_mat_scaled[0, 0] /= alpha
            step_max, eigval_max, nu_max = self.solve_rfo(max_mat_scaled, "max")
            step_max = step_max[0]
            scaled_min_eigenvalues = np.zeros_like(eigvals)
            scaled_min_eigenvalues[:-1] = eigvals[min_indices] / alpha
            min_mat_scaled = min_mat.copy()
            min_mat_scaled[min_diag_indices] = scaled_min_eigenvalues
            step_min, eigval_min, nu_min = self.solve_rfo(min_mat_scaled, "min")

            # As of Eq. (8a) of [4] max_eigval and min_eigval also
            # correspond to:
            # max_eigval = -forces_trans[self.root] * max_step
            # min_eigval = -forces_trans[min_indices].dot(min_step)

            # Create the full PRFO step
            prfo_step = np.zeros_like(forces)
            prfo_step[self.root] = step_max
            prfo_step[min_indices] = step_min
            prfo_norm = np.linalg.norm(prfo_step)
            print("prfo_norm", prfo_norm)

            inside_trust = prfo_norm <= self.trust_radius
            if inside_trust:
                print("PASST")
                break

            # Derivative of the squared step w.r.t. alpha
            # max subspace
            dstep2_dalpha_max = (2*eigval_max/(1+step_max**2 * alpha)
                                 * forces_trans[self.root]**2
                                 / (eigvals[self.root] - eigval_max * alpha)**3
            )
            print("d2da_max", dstep2_dalpha_max)
            # min subspace
            dstep2_dalpha_min = (2*eigval_min/(1+step_min.dot(step_min) * alpha)
                                 * np.sum(forces_trans[min_indices]**2
                                          / (eigvals[min_indices] - eigval_min * alpha)**3
                                 )
            )
            print("d2da_min", dstep2_dalpha_min)
            dstep2_dalpha = dstep2_dalpha_max + dstep2_dalpha_min
            print("d2da", dstep2_dalpha)
            # Update alpha
            alpha_step = (2*(self.trust_radius*prfo_norm - prfo_norm**2)
                          / dstep2_dalpha
            )
            print("alpha_step", alpha_step)
            # alpha += alpha_step
            # alpha *= 0.8#1.1

        predicted_energy_change = 1/2 * eigval_max / nu_max**2 + eigval_min / nu_min**2
        self.predicted_energy_changes.append(predicted_energy_change)
        print("predicted_energy_change", predicted_energy_change)
        # Right now the step is still given in the Hessians eigensystem. We
        # transform it back now.
        step = eigvecs.dot(prfo_step)
        # norm = np.linalg.norm(step)
        # if norm > self.max_size:
            # self.log(f"norm(step, unscaled)={norm:.6f}")
            # self.log("Scaling down step")
            # step = self.max_size * step / norm
            # norm = np.linalg.norm(step)
        # self.log(f"norm(step)={norm:6f}")

        self.log("")
        return step
