# [1] https://aip.scitation.org/doi/pdf/10.1063/1.4905665?class=pdf

import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.intcoords.setup import get_bond_mat
from pysisyphus.optimizers.restrict_step import scale_by_max_step


class StabilizedQNMethod(Optimizer):
    def __init__(
        self,
        geometry,
        alpha=0.5,
        alpha_max=1,
        alpha_stretch=0.5,
        alpha_stretch_max=1,
        eps=1e-4,
        hist_max=10,
        E_thresh=1e-6,
        bio=True,
        trust_radius=0.1,
        linesearch=True,
        **kwargs,
    ):
        super().__init__(geometry, **kwargs)

        self.alpha = alpha
        self.alpha_max = alpha_max
        self.alpha_stretch = alpha_stretch
        self.alpha_stretch_max = alpha_stretch_max
        self.eps = eps
        self.hist_max = int(hist_max)
        self.E_thresh = E_thresh
        self.bio = bio
        self.trust_radius = trust_radius
        self.linesearch = linesearch

        self.log(
            f"Keeping at most information from {self.hist_max} " "previous cycles."
        )

        self.alpha_start = self.alpha
        self.alpha_stretch_start = self.alpha_stretch

        self.gradients_for_precon = list()
        self.coords_for_precon = list()
        self.stretch_proj_signs = list()

    def prepare_opt(self):
        pass

    @property
    def n_hist(self):
        return len(self.steps_normed)

    def bio_mode(self, gradient):
        bond_mat = get_bond_mat(self.geometry)
        atom_num = len(self.geometry.atoms)
        bond_vec_empty = np.zeros((atom_num, 3))
        c3d = self.geometry.coords3d
        bond_vectors = list()
        unique_bonds = set(
            [(min(m, k), max(m, k)) for m, k in zip(*np.where(bond_mat == True))]
        )
        for m, k in unique_bonds:
            bm = bond_vec_empty.copy()
            bm[k] = c3d[k] - c3d[m]
            bm[m] = -bm[k]
            bond_vectors.append(bm)
        # The 3d array is reshaped to 2d, so its last dimension matches
        # the gradient's dimension. LHS and RHS are accessible by simple
        # dot products (see [1] (27)).
        bond_vectors = np.array(bond_vectors).reshape(-1, gradient.size)
        lhs = bond_vectors.dot(gradient)
        self.stretch_proj_signs.append(np.sign(lhs).astype(int))
        # Bond-vector overlap matrix.
        rhs = bond_vectors.dot(bond_vectors.T)
        coeffs = np.linalg.solve(rhs, lhs)
        stretch_gradient = np.einsum("i,ij->j", coeffs, bond_vectors)
        remainder_gradient = gradient - stretch_gradient
        return stretch_gradient, remainder_gradient

    def adjust_alpha_stretch(self):
        try:
            # -1 when signs changed, 1 when sign is constant
            mult = self.stretch_proj_signs[-1] * self.stretch_proj_signs[-2]
        except IndexError:
            self.log("Can't update alpha_stretch yet!")
            return
        except ValueError:
            self.log(
                "Number of bonds found changed between cycles! Skipping "
                "adjustment of alpha_stretch."
            )
            return
        constant_signs = np.sum(mult == 1)
        fraction = constant_signs / mult.size
        if fraction >= (2 / 3):
            self.alpha_stretch = min(self.alpha_stretch_max, 1.1 * self.alpha_stretch)
            msg = "Increased"
        else:
            self.alpha_stretch *= 1 / 1.1
            msg = "Decreased"
        self.log(f"{msg} alpha_stretch to {self.alpha_stretch:.6f}")

    def adjust_alpha(self, gradient, precon_gradient):
        # Overlap between full and preconditioned gradient.
        grad_ovlp = (
            gradient.dot(precon_gradient)
            / np.linalg.norm(gradient)
            / np.linalg.norm(precon_gradient)
        )
        if grad_ovlp > 0.2:
            self.alpha = min(self.alpha_max, 1.1 * self.alpha)
        else:
            self.alpha *= 0.85
        self.log(
            "Overlap between full gradient and preconditioned gradient "
            f"is {grad_ovlp:.4f}. New alpha={self.alpha:.4f}."
        )

    def precondition_gradient(self, gradient, steps, grad_diffs, eps):
        # Construct overlap matrix S between normalized steps
        #
        step_norms = np.linalg.norm(steps, axis=1)
        steps_normed = steps / step_norms[:, None]
        S = np.einsum("ij,jk->ik", steps_normed, steps_normed.T)
        # Determine significant subspace by checking the eigenvalues of
        # the overlap matrix.
        w, v = np.linalg.eigh(S)
        # Significant subspace indices
        significant = (w / w.max()) > eps
        ndim = np.sum(significant)
        self.log(f"ndim = {ndim}")
        # Continue only with eigenvalues and -vectors in the significant
        # subspace.
        eigvals = w[significant]
        eigvecs = v[:, significant]

        # Transform steps and gradient differences to significant subspace
        norm_fac = 1 / eigvals ** (1 / 2)
        # Eq. (8) in [1]
        steps_normed_sub = norm_fac * np.einsum("ki,kj->ji", eigvecs, steps_normed)

        # Eq. (11) in [1]
        eigvecs_weighted = eigvecs / np.array(step_norms)[:, None]
        grad_diffs_sub = norm_fac * np.einsum("ki,kj->ji", eigvecs_weighted, grad_diffs)

        # Assemble approximate hessian
        hess_approx = np.einsum("ij,ik->jk", grad_diffs_sub, steps_normed_sub)
        hess_approx = (hess_approx + hess_approx.T) / 2
        hess_w, hess_v = np.linalg.eigh(hess_approx)

        # Calculate step directions and gradient difference in the original
        # space. Eq. (15)
        proj_v = np.einsum("ki,jk->ji", hess_v, steps_normed_sub)
        proj_dg = np.einsum("ki,jk->ji", hess_v, grad_diffs_sub)

        residuals = np.linalg.norm(proj_dg - hess_w * proj_v, axis=0)
        eigvals_mod = np.sqrt(hess_w ** 2 + residuals ** 2)

        # precon_grad = np.einsum("i,j,ij,ij->i", cur_grad, 1/eigvals_mod, proj_v, proj_v)
        precon_grad = np.einsum(
            "i,j,ij,ij->i", gradient, 1 / eigvals_mod, proj_v, proj_v
        )

        projector = proj_v.dot(proj_v.T)
        perp_projector = np.eye(gradient.size) - projector
        perp_grad = perp_projector.dot(gradient)
        # Alternative formulation
        # perp_grad_ = gradient - (proj_v.T.dot(gradient) * proj_v).sum(axis=1)
        # np.testing.assert_allclose(perp_grad_, perp_grad, atol=1e-16)

        self.adjust_alpha(gradient, precon_grad)
        tot_precon_gradient = precon_grad + self.alpha * perp_grad
        return tot_precon_gradient

    def optimize(self):
        gradient = self.geometry.gradient
        energy = self.geometry.energy
        self.forces.append(-gradient)
        self.energies.append(energy)
        self.log(f"norm(forces)={np.linalg.norm(gradient):.4e}")

        if self.bio:
            stretch_gradient, remainder_gradient = self.bio_mode(gradient)
            self.adjust_alpha_stretch()
            # Steepest descent against the stretch_gradient
            stretch_step = -self.alpha_stretch * stretch_gradient
            new_coords = self.geometry.coords + stretch_step
            self.geometry.coords = new_coords
            # Use only the remaining gradient for the rest of this method
            gradient = remainder_gradient

        self.gradients_for_precon.append(gradient)
        self.coords_for_precon.append(self.geometry.coords.copy())

        if len(self.coords_for_precon) > 2:
            steps = np.diff(self.coords_for_precon, axis=0)[-self.hist_max :]
            grad_diffs = np.diff(self.gradients_for_precon, axis=0)[-self.hist_max :]

            self.log(
                "Preconditioning gradient with information from "
                f"{len(steps)+1} previous cycles."
            )
            precon_grad = self.precondition_gradient(
                gradient, steps, grad_diffs, self.eps
            )
            print("ohoh")
            step = -precon_grad
        else:
            self.log("Took pure steepest descent step.")
            step = self.alpha * -gradient
            scale_by_max_step(step, self.trust_radius)

        if self.linesearch:
            # Calculate energy at new geometry
            new_coords = self.geometry.coords + step
            _ = self.geometry.get_energy_and_forces_at(new_coords)
            new_energy = _["energy"]
            delta_energy = new_energy - energy
            self.log(
                f"Current energy is {energy:.6f} au. New energy is "
                f"{new_energy:.6f} au, ΔE={delta_energy:.6f} au."
            )

            energy_rise_too_big = new_energy > (energy + self.E_thresh)
            alpha_still_big_enough = self.alpha > (self.alpha_start / 10)
            alpha_stre_still_big_enough = self.alpha_stretch > (
                self.alpha_stretch_start / 10
            )
            if (
                energy_rise_too_big
                and alpha_still_big_enough
                and alpha_stre_still_big_enough
            ):
                self.log(f"Energy increased by {delta_energy:.6f} au")
                self.gradients_for_precon = self.gradients_for_precon[-2:-1]
                self.coords_for_precon = self.coords_for_precon[-2:-1]
                self.log("Resetted history.")
                self.alpha /= 2
                self.alpha_stretch /= 2

                self.log(f"Decreased alpha to {self.alpha}")
                self.log("Reverting bad step")
                return None

        return step
