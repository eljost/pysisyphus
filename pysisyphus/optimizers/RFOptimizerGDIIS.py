#!/usr/bin/env python3

# [1] http://aip.scitation.org/doi/10.1063/1.1515483 Optimization review
# [2] https://doi.org/10.1063/1.450914 Trust region method
# [3] 10.1007/978-0-387-40065-5 Numerical optimization
# [4] 10.1007/s00214-016-1847-3 Explorations of some refinements


import numpy as np

from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer
from pysisyphus.optimizers.diis import gdiis


class RFOptimizer(HessianOptimizer):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.error_vecs = list()

    def optimize(self):
        gradient = self.geometry.gradient
        self.forces.append(-self.geometry.gradient)
        self.energies.append(self.geometry.energy)

        if self.cur_cycle > 0:
            self.update_trust_radius()
            self.update_hessian()

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

        self.error_vecs.append(step.copy())
        # if len(self.error_vecs) > 5:
        # # if False:
            # steps_ = self.steps.copy()
            # steps_.append(step)
            # gdiis_kwargs = {
                # "coords": self.coords,
                # "steps": steps_,
                # # "steps": self.steps,
                # "error_vecs": self.error_vecs,
                # "ref_step": step.copy(),
            # }
            # # gdiis_step = gdiis(**gdiis_kwargs)
            # gdiis_coords = gdiis(**gdiis_kwargs)
            # step = gdiis_coords - self.geometry.coords

        step_norm = np.linalg.norm(step)
        self.log(f"norm(step,unscaled)={np.linalg.norm(step):.6f}")
        # Restrict step_norm to the current trust radius
        if step_norm > self.trust_radius:
            self.log(f"step-norm exceeds trust radius; scaling step.")
            step = step / step_norm * self.trust_radius
        self.log(f"norm(step)={np.linalg.norm(step):.6f}")

        predicted_change = step.dot(gradient) + 0.5 * step.dot(self.H).dot(step)
        self.predicted_energy_changes.append(predicted_change)

        return step
