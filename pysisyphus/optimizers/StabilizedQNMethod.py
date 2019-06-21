#!/usr/bin/env python3

# [1] https://aip.scitation.org/doi/pdf/10.1063/1.4905665?class=pdf

import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer


class StabilizedQNMethod(Optimizer):

    def __init__(self, geometry, alpha=0.1, eps=1e-4, **kwargs):
        super().__init__(geometry, **kwargs)

        self.alpha = alpha
        self.eps = eps

        self.step_norms = list()
        self.steps_normed = list()
        self.grad_diffs = list()

    def prepare_opt(self):
        pass

    @property
    def n_hist(self):
        return len(self.steps_normed)

    def get_significant_subspace(self):
        size = len(self.steps_normed)
        # overlap matrix
        S = np.zeros((size, size))
        for i, step_i in enumerate(self.steps_normed):
            for j, step_j in enumerate(self.steps_normed):
                S[i, j] = step_i.dot(step_j)
        w, v = np.linalg.eigh(S)
        # Significant subspace
        significant = (w/w.max()) > self.eps
        ndim = np.sum(significant)
        eigvals = w[significant]
        eigvecs = v[:,significant]
        norm_fac = 1/eigvals**(1/2)
        # Eq. (8) in [1]
        sig_subspace = norm_fac * \
                       np.einsum("ki,kj->ji", eigvecs, self.steps_normed)
        eigvecs_weighted = eigvecs / np.array(self.step_norms)[:, None]
        # Eq. (11) in [1]
        sig_grad_diffs = norm_fac * \
                         np.einsum("ki,kj->ji", eigvecs_weighted, self.grad_diffs)
        hess_approx = np.einsum("ij,ik->jk", sig_grad_diffs, sig_subspace)
        hess_approx = (hess_approx + hess_approx.T) / 2
        hess_w, hess_v = np.linalg.eigh(hess_approx)

        import pdb; pdb.set_trace()
        pass

    def optimize(self):
        forces = self.geometry.forces
        energy = self.geometry.energy
        self.forces.append(forces)
        self.energies.append(energy)

        if len(self.forces) > 1:
            grad_diff = -self.forces[-1] + self.forces[-2]
            self.grad_diffs.append(grad_diff)

        if len(self.steps_normed) > 3:
            self.get_significant_subspace()

        dir_ = forces / np.linalg.norm(forces)
        step = forces
        if np.linalg.norm(step) > 0.1:
            step = self.alpha * dir_


        step_norm = np.linalg.norm(step)
        self.step_norms.append(step_norm)
        step_normed = step / step_norm
        self.steps_normed.append(step_normed)

        return step
        # import pdb; pdb.set_trace()
