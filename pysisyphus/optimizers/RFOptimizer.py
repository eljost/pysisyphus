#!/usr/bin/env python3

# [1] http://aip.scitation.org/doi/10.1063/1.1515483 Optimization review
# [2] https://doi.org/10.1063/1.450914 Trust region method
# [3] 10.1007/978-0-387-40065-5 Numerical optimization
# [4] 10.1007/s00214-016-1847-3 Explorations of some refinements


import numpy as np
from scipy.optimize import minimize

from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.hessian_updates import bfgs_update


class RFOptimizer(Optimizer):

    def __init__(self, geometry, trust_radius=0.3, update_trust=True,
                 **kwargs):
        super().__init__(geometry, **kwargs)

        self.trust_radius = trust_radius
        self.update_trust = update_trust

        self.min_trust_radius = 0.25*self.trust_radius
        self.max_trust_radius = 2
        self.predicted_energy_changes = list()
        self.actual_energy_changes = list()

    def prepare_opt(self):
        self.H = self.geometry.get_initial_hessian()

    def update_trust_radius(self, coeff, last_step_norm):
        # Nocedal, Numerical optimization Chapter 4, Algorithm 4.1
        if coeff < 0.25:
            self.trust_radius = max(self.trust_radius/4,
                                    self.min_trust_radius)
            self.log("Decreasing trust radius.")
        # Only increase trust radius if last step norm was at least 80% of it
        # See [5], Appendix, step size and direction control
        elif coeff > 0.75 and (last_step_norm >= .8*self.trust_radius):
            self.trust_radius = min(self.trust_radius*2,
                                    self.max_trust_radius)
            self.log("Increasing trust radius.")
        else:
            self.log("Keeping current trust radius")
            return
        self.log(f"Updated trust radius: {self.trust_radius:.6f}")

    def optimize(self):
        gradient = self.geometry.gradient
        self.forces.append(-self.geometry.gradient)
        self.energies.append(self.geometry.energy)

        if self.cur_cycle > 0:
            predicted_energy_change = self.predicted_energy_changes[-1]
            actual_energy_change = self.energies[-1] - self.energies[-2]
            if self.update_trust:
                coeff = actual_energy_change / predicted_energy_change
                self.log(f"Predicted change: {predicted_energy_change:.4e} au")
                self.log(f"Actual change: {actual_energy_change:.4e} au")
                self.log(f"Coefficient: {coeff:.2%}")
                last_step_norm = np.linalg.norm(self.steps[-1])
                self.update_trust_radius(coeff, last_step_norm)
            dx = self.steps[-1]
            dg = -(self.forces[-1] - self.forces[-2])
            dH, _ = bfgs_update(self.H, dx, dg)
            self.H += dH

        H = self.H
        if self.geometry.internal:
            H_proj = self.geometry.internal.project_hessian(self.H)
            # Symmetrize hessian, as the projection probably breaks it.
            H = (H_proj + H_proj.T) / 2
        # Eq. (56) in [1]
        H_aug = np.bmat(
                    ((H, gradient[:, None]),
                     (gradient[None, :], [[0]]))
        )
        eigvals, eigvecs = np.linalg.eigh(H_aug)
        # Select eigenvector corresponding to smallest eigenvalue. Eigen-
        # values and -vectors are sorted, so we take the first one.
        aug_step = eigvecs[:, 0]
        # aug_step is currently a matrix. Convert it to an array.
        aug_step = np.asarray(aug_step).flatten()
        lambda_ = aug_step[-1]
        # Scale aug_step so the last element equals 1, then neglect the last
        # item.
        step = aug_step[:-1] / lambda_

        step_norm = np.linalg.norm(step)
        # Restrict step_norm to the current trust radius
        if step_norm > self.trust_radius:
            step = step / step_norm * self.trust_radius
        self.log(f"norm(step): {np.linalg.norm(step):.4f}")

        predicted = step.dot(gradient) + 0.5 * step.dot(H).dot(step)
        self.predicted_energy_changes.append(predicted)

        return step
