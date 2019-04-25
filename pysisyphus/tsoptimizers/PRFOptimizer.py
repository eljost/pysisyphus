#!/usr/bin/env python3

# See [1] https://pubs.acs.org/doi/pdf/10.1021/j100247a015
#         Banerjee, 1985
#     [2] https://aip.scitation.org/doi/abs/10.1063/1.2104507
#         Heyden, 2005
#     [3] https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.540070402
#         Baker, 1985
# TODO: 10.1007/s002140050387


import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer


class PRFOptimizer(Optimizer):
    """Optimizer to find first-order saddle points."""

    def __init__(self, geometry, root=0, max_step_length=.3, recalc_hess=None,
                 **kwargs):
        super().__init__(geometry, **kwargs)

        self.root = int(root)
        self.max_step_length = max_step_length
        self.recalc_hess = recalc_hess

        self.H = None
        self.ts_mode = None

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
            # dx = self.coords[-1] - self.coords[-2]
            dx = self.steps[-1]
            self.H += self.bofill_update(self.H, dx, dg)
            self.log("Did Bofill hessian update.")

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
        min_mat = np.bmat((
            (np.diag(eigvals[min_indices]), -forces_trans[min_indices,None]),
            (-forces_trans[None,min_indices], [[0]])
        ))

        # Scale eigenvectors of the largest (smallest) eigenvector
        # of max_mat (min_mat) so the last item equals to 1.
        # An alternative route is given in [3] Sec. II.3 by solving
        # a quadratic equation for lambda_max.
        max_eigvals, max_evecs = np.linalg.eigh(max_mat)
        # Eigenvalues and -values are sorted, so we just use the last
        # eigenvector corresponding to the biggest eigenvalue.
        max_step = max_evecs.T[-1]
        lambda_max = max_step[-1]
        self.log(f"lambda_max={lambda_max:.4e}")
        # Drop last element
        max_step = max_step[:-1] / lambda_max
        max_shift = max_eigvals[-1] / lambda_max
        self.log(f"shift_maximize={max_shift:.4e}")

        min_eigvals, min_evecs = np.linalg.eigh(min_mat)
        # Again, as everything is sorted we use the (smalelst) first eigenvalue.
        min_step = np.asarray(min_evecs.T[0]).flatten()
        lambda_min = min_step[-1]
        self.log(f"lambda_min={lambda_min:.4e}")
        # Drop last element
        min_step = min_step[:-1] / lambda_min
        min_shift = min_eigvals[0] / lambda_min
        self.log(f"shift_minimize={min_shift:.4e}")

        # Create the full PRFO step
        prfo_step = np.zeros_like(forces)
        prfo_step[self.root] = max_step[0]
        prfo_step[min_indices] = min_step
        # Right now the step is still given in the Hessians eigensystem. We
        # transform it back now.
        step = eigvecs.dot(prfo_step)
        norm = np.linalg.norm(step)
        prfo_dir = step / norm
        if norm > self.max_step_length:
            self.log(f"norm(step, unscaled)={norm:.6f}")
            self.log("Scaling down step")
            step = self.max_step_length * step / norm
            norm = np.linalg.norm(step)
        self.log(f"norm(step)={norm:6f}")

        self.log("")
        return step
