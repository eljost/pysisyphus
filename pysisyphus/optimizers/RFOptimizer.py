#!/usr/bin/env python3

# [1] http://aip.scitation.org/doi/10.1063/1.1515483 Optimization review
# [2] https://doi.org/10.1063/1.450914 Trust region method
# [3] 10.1007/978-0-387-40065-5 Numerical optimization
# [4] 10.1007/s00214-016-1847-3 Explorations of some refinements


import numpy as np

from pysisyphus.helpers import rms
from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer
from pysisyphus.optimizers.gdiis import gdiis, gediis


class RFOptimizer(HessianOptimizer):

    def __init__(self, geom, line_search=True, gediis=False, gdiis=True,
                 gdiis_thresh=2.5e-3, gediis_thresh=1e-2,
                 *args, **kwargs):
        super().__init__(geom, *args, **kwargs)

        self.line_search = line_search
        self.gediis = gediis
        self.gdiis = gdiis
        self.gdiis_thresh = gdiis_thresh  # Will be compared to rms(step)
        self.gediis_thresh = gediis_thresh  # Will be compared to rms(forces)

    def optimize(self):
        energy, gradient, H, big_eigvals, big_eigvecs = self.housekeeping()

        org_grad = gradient.copy()

        dim_ = big_eigvals.size + 1
        def get_step(gradient, eigvals, eigvecs):
            gradient_ = big_eigvecs.T @ gradient
            H_aug = np.zeros((dim_, dim_))
            H_aug[:dim_-1,:dim_-1] = np.diag(big_eigvals)
            H_aug[-1,:-1] = gradient_
            H_aug[:-1,-1] = gradient_
            step_, eigval, nu = self.solve_rfo(H_aug, "min")
            # Transform back to original basis
            step = big_eigvecs @ step_
            return step

        ref_step = get_step(gradient, big_eigvals, big_eigvecs)
        # Right now we have everything in place to check for convergence.
        # If all values are below the thresholds there is no need to do additional
        # inter/extrapolations.
        if self.check_convergence(ref_step):
            self.log("Convergence achieved! Skipping inter/extrapolation.")
            return  ref_step
        step = ref_step

        rms_forces = rms(self.forces[-1])
        rms_step = rms(ref_step)
        can_diis = rms_step <= self.gdiis_thresh
        can_gediis = rms_forces <= self.gediis_thresh
        diis_result = None
        ip_gradient = None
        if self.gdiis and can_diis:
            err_vecs = -np.array(self.forces)
            diis_result = gdiis(err_vecs, self.coords, self.forces, ref_step)
        # Don't try GEDIIS if GDIIS failed
        elif self.gediis and can_gediis:
        # Try GEDIIS if GDIIS failed
        # if self.gediis and can_gediis and (diis_result == None):
            diis_result = gediis(self.coords, self.energies, self.forces, hessian=H)

        try:
            ip_coords = diis_result.coords
            tmp_geom = self.geometry.copy()
            tmp_geom.coords = diis_result.coords
            ip_step = tmp_geom - self.geometry
            ip_gradient = -diis_result.forces.copy()
        # When diis_result is None
        except AttributeError:
            self.log("GDIIS didn't succeed.")
        except ValueError:
            # This will be raised if the tmp_geom will have a different
            # number of coordinates because some (primitives) couldn't
            # be defined.
            diis_result = None

        can_linesearch = (diis_result is None) and self.line_search and (self.cur_cycle > 0)
        if can_linesearch:
            ip_energy, ip_gradient, ip_coords, ip_step = self.poly_line_search()

        if ip_gradient is not None:
            # Project interpolated gradient if necessary
            if self.geometry.coord_type == "redund":
                ip_gradient = self.geometry.internal.project_vector(ip_gradient)

            self.geometry.coords = ip_coords
            self.forces[-1] = -ip_gradient
            try:
                self.energies[-1] = ip_energy
            # GDIIS doesn't produce an energy
            except UnboundLocalError:
                self.log("Skipping energy update after fitting.")
            self.coords[-1] = ip_coords.copy()
            self.cart_coords[-1] = self.geometry.cart_coords.copy()
            self.steps[-1] = ip_step

            source = diis_result.type if diis_result else "line search"
            self.log(f"Got inter/extra-polated geometry from {source}")
            self.log(f"Calculating second step from inter/extra-polated gradient.")
            if self.geometry.coord_type == "redund":
                ip_gradient = self.geometry.internal.project_vector(ip_gradient)
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
