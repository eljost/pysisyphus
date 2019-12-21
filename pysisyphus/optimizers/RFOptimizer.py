#!/usr/bin/env python3

# [1] http://aip.scitation.org/doi/10.1063/1.1515483 Optimization review
# [2] https://doi.org/10.1063/1.450914 Trust region method
# [3] 10.1007/978-0-387-40065-5 Numerical optimization
# [4] 10.1007/s00214-016-1847-3 Explorations of some refinements


import numpy as np

from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer
from pysisyphus.optimizers.gdiis import gdiis, gediis
from pysisyphus.optimizers.interpolate_extrapolate import interpolate_extrapolate


class RFOptimizer(HessianOptimizer):

    def __init__(self, geom, line_search=True, gediis=False, gdiis=True,
                 *args, **kwargs):
        super().__init__(geom, *args, **kwargs)

        self.line_search = line_search
        self.gediis = gediis
        self.gdiis = gdiis

    def optimize(self):
        gradient = self.geometry.gradient
        self.forces.append(-self.geometry.gradient)
        self.energies.append(self.geometry.energy)

        org_grad = gradient.copy()
        if self.cur_cycle > 0:
            # TODO: skip hessian update in the beginning when the gradients
            # are big?
            self.update_trust_radius()
            self.update_hessian()

        H = self.H
        if self.geometry.internal:
            # Shift eigenvalues of orthogonal part to high values, so they
            # don't contribute to the actual step.
            H_proj = self.geometry.internal.project_hessian(H)
            # Symmetrize hessian, as the projection probably breaks it?!
            H = (H_proj + H_proj.T) / 2

        eigvals, eigvecs = np.linalg.eigh(H)

        # Calculate step in basis of eigenvectors of the hessian.
        big_eigvals, big_eigvecs = self.filter_small_eigvals(eigvals, eigvecs)
        dim_ = big_eigvals.size + 1
        def get_step(gradient, eigvals, eigvecs):
            gradient_ = big_eigvecs.T @ gradient
            # H_aug = np.zeros((dim_, dim_))
            # H_aug[:dim_-1,:dim_-1] = np.diag(big_eigvals)
            # H_aug[-1,:-1] = gradient_
            # H_aug[:-1,-1] = gradient_
            H_aug = np.array(
                np.bmat((
                    (np.diag(big_eigvals), gradient_[:, None]),
                    (gradient_[None,:], [[0]])
                ))
            )
            step_, eigval, nu = self.solve_rfo(H_aug, "min")
            # Transform back to original basis
            step = big_eigvecs @ step_
            return step

        ref_step = get_step(gradient, big_eigvals, big_eigvecs)
        step = ref_step

        gediis_thresh = 1e-2
        gdiis_thresh = 2.5e-3
        can_gediis = np.sqrt(np.mean(self.forces[-1]**2)) < gediis_thresh
        can_diis = np.sqrt(np.mean(ref_step**2)) < gdiis_thresh
        diis_result = None
        ip_gradient = None
        if self.gdiis and can_diis:
            err_vecs = -np.array(self.forces)
            diis_result = gdiis(err_vecs, self.coords, self.forces, ref_step)
        elif self.gediis and can_gediis:
            diis_result = gediis(self.coords, self.energies, self.forces)

        if diis_result:
            tmp_geom = self.geometry.copy()
            try:
                tmp_geom.coords = diis_result.coords
                diis_step = tmp_geom - self.geometry
                self.geometry.coords = diis_result.coords
                self.forces[-1] = diis_result.forces
                # self.energies[-1] = y
                self.coords[-1] = self.geometry.coords.copy()
                self.cart_coords[-1] = self.geometry.cart_coords.copy()
                self.steps[-1] = diis_step
                ip_gradient = -diis_result.forces
            except ValueError:
                # This will be raised if the tmp_geom will have a different
                # number of coordinates because some (primitives) couldn't
                # be defined.
                diis_result = None

        can_linesearch = (diis_result is None) and self.line_search and (self.cur_cycle > 0)
        if can_linesearch:
            ip_gradient = self.poly_line_search()

        if ip_gradient is not None:
            step = get_step(ip_gradient, big_eigvals, big_eigvecs)

        step_norm = np.linalg.norm(step)
        self.log(f"norm(step,unscaled)={np.linalg.norm(step):.6f}")
        # Restrict step_norm to the current trust radius
        if step_norm > self.trust_radius:
            self.log(f"step-norm exceeds trust radius; scaling step.")
            step = step / step_norm * self.trust_radius
        self.log(f"norm(step)={np.linalg.norm(step):.6f}")

        # quadratic_prediction = step @ gradient + 0.5 * step @ self.H @ step
        quadratic_prediction = step @ org_grad + 0.5 * step @ self.H @ step
        rfo_prediction = quadratic_prediction / (1 + step @ step)
        self.predicted_energy_changes.append(rfo_prediction)

        return step
