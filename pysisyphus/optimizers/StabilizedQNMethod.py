#!/usr/bin/env python3

# [1] https://aip.scitation.org/doi/pdf/10.1063/1.4905665?class=pdf

import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.guess_hessians import get_bond_mat


class StabilizedQNMethod(Optimizer):

    def __init__(self, geometry, alpha=0.1, alpha_max=1,
                 eps=1e-4, hist_max=10, E_thresh=1e-6, bio=True, **kwargs):
        super().__init__(geometry, **kwargs)

        self.alpha = alpha
        self.alpha_max = alpha_max
        self.eps = eps
        self.hist_max = int(hist_max)
        self.E_thresh = E_thresh
        self.bio = bio

        self.log(f"Keeping at most information from {self.hist_max} "
                  "previous cycles.")

        self.alpha_start = self.alpha

        self.step_norms = list()
        self.steps_normed = list()
        self.grad_diffs = list()

    def prepare_opt(self):
        pass

    @property
    def n_hist(self):
        return len(self.steps_normed)

    def bio_mode(self):
        bond_mat = get_bond_mat(self.geometry)
        from scipy.spatial.distance import squareform
        sq = squareform(bond_mat)
        zero = np.zeros(3)
        atom_num = len(self.geometry.atoms)
        bond_vec_empty = np.zeros((atom_num, 3))
        c3d = self.geometry.coords3d
        bond_vectors = list()
        # unique_bonds = set([(m, k) for m, k in zip(*np.where(bond_mat == True))])
        import pdb; pdb.set_trace()
        for m, k in zip(*np.where(bond_mat == True)):
            print("bond", m, k)
            bm = bond_vec_empty.copy()
            bm[k] = c3d[k] - c3d[m]
            bm[m] = -bm[k]
            bond_vectors.append(bm)
        bond_vectors = np.array(bond_vectors)
        import pdb; pdb.set_trace()
        pass

    def precondition_gradient(self):
        # Construct overlap matrix S between normalized steps
        #
        # As self.steps_normed is a list and not an array we can't
        # easily transpose it. So we swap 'jk' to 'kj' in the second term
        # of the einsum call.
        # The call with transposed self.steps_normed would be
        #    np.einsum("ij,jk->ik", self.steps_normed, self.steps_normed.T)
        S = np.einsum("ij,kj->ik", self.steps_normed, self.steps_normed)
        # Determine significant subspace by checking the eigenvalues of
        # the overlap matrix.
        w, v = np.linalg.eigh(S)
        # Significant subspace indices
        significant = (w/w.max()) > self.eps
        ndim = np.sum(significant)
        self.log(f"ndim = {ndim}")
        # Continue only with eigenvalues and -vectors in the significant
        # subspace.
        eigvals = w[significant]
        eigvecs = v[:,significant]

        # Transform steps and gradient differences to significant subspace
        norm_fac = 1/eigvals**(1/2)
        # Eq. (8) in [1]
        steps_normed_sub = norm_fac * \
                           np.einsum("ki,kj->ji", eigvecs, self.steps_normed)

        # Eq. (11) in [1]
        eigvecs_weighted = eigvecs / np.array(self.step_norms)[:, None]
        grad_diffs_sub = norm_fac * \
                         np.einsum("ki,kj->ji", eigvecs_weighted, self.grad_diffs)

        # Assemble approximate hessian
        hess_approx = np.einsum("ij,ik->jk", grad_diffs_sub, steps_normed_sub)
        hess_approx = (hess_approx + hess_approx.T) / 2
        hess_w, hess_v = np.linalg.eigh(hess_approx)

        # Calculate step directions and gradient difference in the original
        # space. Eq. (15)
        proj_v = np.einsum("ki,jk->ji", hess_v, steps_normed_sub)
        proj_dg = np.einsum("ki,jk->ji", hess_v, grad_diffs_sub)

        residuals = np.linalg.norm(proj_dg - hess_w*proj_v, axis=0)
        eigvals_mod = np.sqrt(hess_w**2 + residuals**2)

        cur_grad = -self.forces[-1]
        precon_grad = np.einsum("i,j,ij,ij->i", cur_grad, 1/eigvals_mod, proj_v, proj_v)

        # projector = proj_v @ proj_v.T
        projector = proj_v.dot(proj_v.T)
        perp_projector = np.eye(cur_grad.size) - projector
        perp_grad = perp_projector.dot(cur_grad)
        # Overlap between full and preconditioned gradient. Need to scale
        # alpha appropriately
        grad_ovlp = (cur_grad.dot(precon_grad)
                     / np.linalg.norm(cur_grad) / np.linalg.norm(precon_grad)
        )
        if grad_ovlp > 0.2:
            self.alpha = min(self.alpha_max, 1.1*self.alpha)
        else:
            self.alpha *= 0.85
        self.log( "Overlap between full gradient and preconditioned gradient "
                 f"is {grad_ovlp:.4f}. New alpha={self.alpha:.4f}.")
        tot_precon_gradient = precon_grad + self.alpha*perp_grad
        return tot_precon_gradient

    def optimize(self):
        forces = self.geometry.forces
        energy = self.geometry.energy
        self.forces.append(forces)
        self.energies.append(energy)

        if len(self.forces) > 1:
            grad_diff = -self.forces[-1] + self.forces[-2]
            self.grad_diffs.append(grad_diff)

        # if self.cur_cycle == 2:
            # import pdb; pdb.set_trace()

        if len(self.steps_normed) > 2:
            assert (len(self.steps_normed) == len(self.step_norms)
                    == len(self.grad_diffs))

            self.log( "Preconditioning gradient with information from "
                     f"{len(self.grad_diffs)} previous cycles.")
            precon_grad = self.precondition_gradient()
            step = -precon_grad
        else:
            dir_ = forces / np.linalg.norm(forces)
            step = self.alpha * forces
            if np.linalg.norm(step) > 0.1:
                step = self.alpha * dir_

        new_coords = self.coords[-1] + step
        self.geometry.coords = new_coords
        new_energy = self.geometry.energy
        delta_energy = new_energy - energy
        self.log(f"Current energy is {energy:.6f} au. New energy is "
                 f"{new_energy:.6f} au ΔE={delta_energy:.6f} au.")

        energy_rise_too_big = new_energy > (energy + self.E_thresh)
        alpha_still_big_enough = self.alpha > (self.alpha_start / 10)
        if energy_rise_too_big and alpha_still_big_enough:
            print("ohoh")
            # # self.step_norms = self.step_norms[-1:]
            # # self.steps_normed = self.steps_normed[-1:]
            # # self.grad_diffs = self.grad_diffs[-1:]
            # self.geometry.coords = self.coords[-1].copy()
            # self.geometry.forces = self.forces[-1]
            # self.geometry.energy = self.energies[-1]
            # self.alpha /= 2
            # self.log(f"Decreased alpha to {self.alpha}")
            # self.log("Reverting bad step")
            # return None

        step_norm = np.linalg.norm(step)
        self.step_norms.append(step_norm)
        step_normed = step / step_norm
        self.steps_normed.append(step_normed)

        self.step_norms = self.step_norms[-self.hist_max:]
        self.steps_normed = self.steps_normed[-self.hist_max:]
        self.grad_diffs = self.grad_diffs[-self.hist_max+1:]

        return step
