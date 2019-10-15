#!/usr/bin/env python3

# [1] http://aip.scitation.org/doi/10.1063/1.1515483 Optimization review
# [2] https://doi.org/10.1063/1.450914 Trust region method
# [3] 10.1007/978-0-387-40065-5 Numerical optimization
# [4] 10.1007/s00214-016-1847-3 Explorations of some refinements


import numpy as np

from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer
# from pysisyphus.optimizers.gdiis import gdiis, gediis
from pysisyphus.optimizers.interpolate_extrapolate import interpolate_extrapolate


class RFOptimizer(HessianOptimizer):

    def __init__(self, geom, eins=False, zwei=False, *args, **kwargs):
        super().__init__(geom, *args, **kwargs)

        self.eins = eins
        self.zwei = zwei

    def optimize(self):
        gradient = self.geometry.gradient
        self.forces.append(-self.geometry.gradient)
        self.energies.append(self.geometry.energy)

        if self.cur_cycle > 0:
            # TODO: skip hessian update in the beginning when the gradients
            # are big?
            self.update_trust_radius()
            self.update_hessian()
            if self.eins:
                # print("eins")
                gradient = self.poly_line_search()
                # print(f"\tnorm(grad)={np.linalg.norm(gradient):.6f}")

        H = self.H

        eigvals, _ = np.linalg.eigh(H)
        neg_eigval_inds = eigvals < -self.small_eigval_thresh
        neg_num = neg_eigval_inds.sum()
        neg_eigval_str = np.array2string(eigvals[neg_eigval_inds], precision=6)
        self.log(f"Found {neg_num} negative eigenvalue(s): {neg_eigval_str}")

        small_eigval_inds = (0 > eigvals) & (eigvals < self.small_eigval_thresh)
        small_num = small_eigval_inds.sum()
        small_eigval_str = np.array2string(eigvals[small_eigval_inds], precision=6)
        self.log(f"Found {small_num} small eigenvalue(s): {small_eigval_str}")

        if self.geometry.internal:
            H_proj = self.geometry.internal.project_hessian(self.H)
            # Symmetrize hessian, as the projection probably breaks it?!
            H = (H_proj + H_proj.T) / 2
        elif self.geometry.coord_type == "cart" and (not self.is_cos):
            pass

        # TODO: Neglect small eigenvalues with cartesian coordinates
        # or use eckard projection.

        # Eq. (56) in [1]
        H_aug = np.array(np.bmat(
                            ((H, gradient[:, None]),
                             (gradient[None, :], [[0]]))
        ))
        step, eigval, nu = self.solve_rfo(H_aug, "min")

        interpol_forces = None
        if self.zwei and self.cur_cycle > 0:
            # print("zwei")
            ie_kwargs = {
                # "geom": self.geometry,
                "coords": self.coords,
                "energies": self.energies,
                "forces": self.forces,
                "steps": self.steps,
                "ref_step": step,
                "err_vecs": self.forces,
                "max_vecs": 10,
                # "gediis_thresh": 1e-2,
                # "gdiis_thresh": 2.5e-3,
                "gediis_thresh": -1,
                "gdiis_thresh": -1,
            }
            interpol_step, interpol_forces, interpol_en = interpolate_extrapolate(**ie_kwargs)
        if interpol_forces is not None:
            # print("drei")
            # fit_coords = self.coords[-2] + interpol_step
            # self.geometry.coords = fit_coords
            # self.forces[-1] = interpol_forces
            # self.energies[-1] = interpol_en
            # self.coords[-1] = fit_coords.copy()
            # self.cart_coords[-1] = self.geometry.cart_coords.copy()
            # self.steps[-1] = interpol_step
            # print(f"\tnorm(grad)={np.linalg.norm(interpol_forces):.6f}")

            # fit_coords = self.coords[-2] + interpol_step
            # fit_coords = self.coords[-1] + interpol_step

            # TODO: update step and other saved entries?!
            # self.geometry.coords = fit_coords
            # self.forces[-1] = interpol_forces
            # self.energies[-1] = interpol_en
            # self.coords[-1] = fit_coords.copy()
            # # self.cart_coords[-1] = self.geometry.cart_coords.copy()
            # self.steps[-1] = interpol_step

            H_aug = np.array(np.bmat(
                                ((H, -interpol_forces[:, None]),
                                 (-interpol_forces[None, :], [[0]]))
            ))
            step, eigval, nu = self.solve_rfo(H_aug, "min")
            step = interpol_step + step
            # step = interpol_step + step
        # else:
            # print(f"\tnorm(grad)={np.linalg.norm(gradient):.6f}")

        # can_gediis = np.sqrt(np.mean(self.forces[-1]**2)) < 1e-2
        # try:
            # can_diis = np.sqrt(np.mean(self.steps[-1]**2)) < 2.5e-3
        # except IndexError:
            # can_diis = False
        # # GDIIS check
        # can_gediis = False
        # if self.hybrid and (self.cur_cycle) > 0 and can_diis:
            # # rms_step = np.sqrt(np.mean(self.steps[-1]**2))
            # # print(f"\t\t\trms(step): {rms_step:.6f}")
            # gdiis_kwargs = {
                # "coords": self.coords,
                # "forces": self.forces,
                # "ref_step": step,
            # }
            # diis_result = gdiis(self.forces, **gdiis_kwargs)
        # # GEDIIS check
        # elif self.hybrid and (self.cur_cycle > 0) and can_gediis:
            # diis_result = gediis(self.coords, self.energies, self.forces)
        # else:
            # diis_result = None

        # if diis_result:
            # # Inter-/extrapolate coordinates and forces
            # forces = diis_result.forces
            # energy = self.energies[-1]
            # # print(f"\torg.-energy={energy:.6f}")
            # import pdb; pdb.set_trace()
            # self.geometry.coords = diis_result.coords
            # diis_energy = self.geometry.energy
            # # print(f"\tdiis-energy={diis_energy:.6f}")
            # # print(f"\tdiis-energy is lower? {diis_energy < energy}")
            # # Get new step from DIIS coordinates & forces
            # H_aug = np.array(np.bmat(
                                # ((H, -forces[:, None]),
                                 # (forces[None, :], [[0]]))
            # ))
            # step, eigval, nu = self.solve_rfo(H_aug, "min")
        # # elif self.line_search:
            # # gradient = self.poly_line_search()

        step_norm = np.linalg.norm(step)
        self.log(f"norm(step,unscaled)={np.linalg.norm(step):.6f}")
        # Restrict step_norm to the current trust radius
        if step_norm > self.trust_radius:
            self.log(f"step-norm exceeds trust radius; scaling step.")
            step = step / step_norm * self.trust_radius
        self.log(f"norm(step)={np.linalg.norm(step):.6f}")

        quadratic_prediction = step @ gradient + 0.5 * step @ self.H @ step
        rfo_prediction = quadratic_prediction / (1 + step @ step)
        self.predicted_energy_changes.append(rfo_prediction)

        return step

    # def optimize(self):
        # gradient = self.geometry.gradient
        # self.forces.append(-self.geometry.gradient)
        # self.energies.append(self.geometry.energy)

        # if self.cur_cycle > 0:
            # # TODO: skip hessian update in the beginning when the gradients
            # # are big?
            # self.update_trust_radius()
            # self.update_hessian()
            # if self.line_search:
                # gradient = self.poly_line_search()

        # H = self.H

        # eigvals, eigvectors = np.linalg.eigh(H)
        # neg_eigval_inds = eigvals < -self.small_eigval_thresh
        # neg_num = neg_eigval_inds.sum()
        # neg_eigval_str = np.array2string(eigvals[neg_eigval_inds], precision=6)
        # self.log(f"Found {neg_num} negative eigenvalue(s): {neg_eigval_str}")

        # small_eigval_inds = (0 > eigvals) & (eigvals < self.small_eigval_thresh)
        # small_num = small_eigval_inds.sum()
        # small_eigval_str = np.array2string(eigvals[small_eigval_inds], precision=6)
        # self.log(f"Found {small_num} small eigenvalue(s): {small_eigval_str}")

        # if self.geometry.internal:
            # H_proj = self.geometry.internal.project_hessian(self.H)
            # # Symmetrize hessian, as the projection probably breaks it?!
            # H = (H_proj + H_proj.T) / 2
            # eigvals, eigvectors = np.linalg.eigh(H)
        # elif self.geometry.coord_type == "cart" and (not self.is_cos):
            # eigvals, eigvectors = self.filter_small_eigvals(eigvals, eigvectors)

        # # Transform gradient to the hessian eigenspace
        # gradient_ = eigvectors.T @ gradient

        # # Construct augmented hessian
        # # Eq. (56) in [1]
        # H_aug = np.array(np.bmat(
                            # ((np.diag(eigvals), gradient_[:, None]),
                             # (gradient_[None, :], [[0]]))
        # ))
        # step_, eigval, nu = self.solve_rfo(H_aug, "min")

        # # Transform back to original space
        # step = eigvectors @ step_


        # can_gediis = np.sqrt(np.mean(self.forces[-1]**2)) < 1e-2
        # try:
            # can_diis = np.sqrt(np.mean(self.steps[-1]**2)) < 2.5e-3
        # except IndexError:
            # can_diis = False
        # # GDIIS check
        # can_gediis = True
        # if self.hybrid and (self.cur_cycle) > 0 and can_diis:
            # # rms_step = np.sqrt(np.mean(self.steps[-1]**2))
            # # print(f"\t\t\trms(step): {rms_step:.6f}")
            # gdiis_kwargs = {
                # "coords": self.coords,
                # "forces": self.forces,
                # "ref_step": step,
            # }
            # diis_result = gdiis(self.forces, **gdiis_kwargs)
        # # GEDIIS check
        # elif self.hybrid and (self.cur_cycle > 0) and can_gediis:
            # diis_result = gediis(self.coords, self.energies, self.forces)
        # else:
            # diis_result = None

        # if diis_result:
            # # Inter-/extrapolate coordinates and forces
            # forces = diis_result.forces
            # energy = self.energies[-1]
            # # print(f"\torg.-energy={energy:.6f}")
            # self.geometry.coords = diis_result.coords
            # diis_energy = self.geometry.energy
            # # print(f"\tdiis-energy={diis_energy:.6f}")
            # # print(f"\tdiis-energy is lower? {diis_energy < energy}")
            # # Get new step from DIIS coordinates & forces
            # H_aug = np.array(np.bmat(
                                # ((H, -forces[:, None]),
                                 # (forces[None, :], [[0]]))
            # ))
            # step, eigval, nu = self.solve_rfo(H_aug, "min")

        # step_norm = np.linalg.norm(step)
        # self.log(f"norm(step,unscaled)={np.linalg.norm(step):.6f}")
        # # Restrict step_norm to the current trust radius
        # if step_norm > self.trust_radius:
            # self.log(f"step-norm exceeds trust radius; scaling step.")
            # step = step / step_norm * self.trust_radius
        # self.log(f"norm(step)={np.linalg.norm(step):.6f}")

        # quadratic_prediction = step @ gradient + 0.5 * step @ self.H @ step
        # rfo_prediction = quadratic_prediction / (1 + step @ step)
        # self.predicted_energy_changes.append(rfo_prediction)

        # return step
