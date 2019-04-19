#!/usr/bin/env python3

# See [1] https://pubs.acs.org/doi/pdf/10.1021/j100247a015
#         Banerjee, 1985
#     [2] https://aip.scitation.org/doi/abs/10.1063/1.2104507
#         Heyden, 2005

import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer


class PRFOptimizer(Optimizer):

    def __init__(self, geometry, root=-1, max_size=.1, **kwargs):
        super().__init__(geometry, **kwargs)

        self.root = root
        self.max_size = max_size

    def prepare_opt(self):
        self.H = self.geometry.hessian

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

        if len(self.coords) > 1:
            dg = -(self.forces[-1] - self.forces[-2])
            dx = self.coords[-1] - self.coords[-2]
            self.H += self.bofill_update(self.H, dx, dg)

        eigvals, eigvecs = np.linalg.eigh(self.H)
        neg_eigvals = eigvals < -1e-8
        assert neg_eigvals.sum() >= 1, \
            "Need at least 1 negative eigenvalue for TS optimization."
        eigval_str = np.array2string(eigvals[neg_eigvals], precision=6)
        self.log(f"Found {neg_eigvals.sum()} negative eigenvalue(s): {eigval_str}")

        # Transform to eigensystem of hessian
        forces_trans = eigvecs.T.dot(forces)

        # Maximize energy along first root
        mu = 0
        max_mat = np.array(((eigvals[mu], -forces_trans[mu]),
                           (-forces_trans[mu], 0)))
        # Minimize energy along the other roots
        min_mat = np.bmat((
            (np.diag(eigvals[1:]), -forces_trans[1:,None]),
            (-forces_trans[None,1:], [[0]])
        ))

        # Scale eigenvectors of the largest (smallest) eigenvector
        # of max_mat (min_mat) so the last item is 1.
        max_evals, max_evecs = np.linalg.eigh(max_mat)
        # Eigenvalues and -values are sorted, so we just use the last
        # eigenvector corresponding to the biggest eigenvalue.
        max_step = max_evecs.T[-1]
        lambda_max = max_step[-1]
        max_step = max_step[:-1] / lambda_max

        min_evals, min_evecs = np.linalg.eigh(min_mat)
        # Again, as everything is sorted we use the (smalelst) first eigenvalue.
        min_step = np.asarray(min_evecs.T[0]).flatten()
        lambda_min = min_step[-1]
        min_step = min_step[:-1] / lambda_min

        # Create the full PRFO step
        prfo_step = np.zeros_like(forces)
        prfo_step[0] = max_step[0]
        prfo_step[1:] = min_step
        # Right now the step is still given in the Hessians eigensystem. We
        # transform it back now.
        step = eigvecs.dot(prfo_step)
        norm = np.linalg.norm(step)
        if norm > self.max_size:
            step = self.max_size * step / norm
        return step
