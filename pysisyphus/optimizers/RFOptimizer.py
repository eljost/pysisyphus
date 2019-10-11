#!/usr/bin/env python3

# [1] http://aip.scitation.org/doi/10.1063/1.1515483 Optimization review
# [2] https://doi.org/10.1063/1.450914 Trust region method
# [3] 10.1007/978-0-387-40065-5 Numerical optimization
# [4] 10.1007/s00214-016-1847-3 Explorations of some refinements


import numpy as np

from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer
from pysisyphus.optimizers.gdiis import gdiis


class RFOptimizer(HessianOptimizer):

    def optimize(self):
        gradient = self.geometry.gradient
        self.forces.append(-self.geometry.gradient)
        self.energies.append(self.geometry.energy)

        if self.cur_cycle > 0:
            # TODO: skip hessian update in the beginning when the gradients
            # are big?
            self.update_trust_radius()
            self.update_hessian()
            if self.line_search:
                gradient = self.poly_line_search()

        H = self.H

        eigvals, _ = np.linalg.eigh(H)
        neg_eigval_inds = eigvals < -self.small_eigval_thresh
        neg_num = neg_eigval_inds.sum()
        eigval_str = np.array2string(eigvals[neg_eigval_inds], precision=6)
        self.log(f"Found {neg_num} negative eigenvalue(s): {eigval_str}")

        if self.geometry.internal:
            H_proj = self.geometry.internal.project_hessian(self.H)
            # Symmetrize hessian, as the projection probably breaks it?!
            H = (H_proj + H_proj.T) / 2

        # TODO: Neglect small eigenvalues with cartesian coordinates
        # or use eckard projection.

        # Eq. (56) in [1]
        H_aug = np.array(np.bmat(
                            ((H, gradient[:, None]),
                             (gradient[None, :], [[0]]))
        ))
        step, eigval, nu = self.solve_rfo(H_aug, "min")

        # if self.cur_cycle > 0:
            # gdiis_kwargs = {
                # "coords": self.coords,
                # "forces": self.forces,
                # "ref_step": step,
            # }
            # gdiis_result = gdiis(self.forces, **gdiis_kwargs)
            # if gdiis_result:
                # # Inter-/extrapolate coordinates and forces
                # forces = gdiis_result.forces
                # self.geometry.coords = gdiis_result.coords
                # # Get new step from DIIS coordinates & forces
                # H_aug = np.array(np.bmat(
                                    # ((H, -forces[:, None]),
                                     # (forces[None, :], [[0]]))
                # ))
                # step, eigval, nu = self.solve_rfo(H_aug, "min")

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
