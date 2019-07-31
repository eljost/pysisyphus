#!/usr/bin/env python3

# [1] http://aip.scitation.org/doi/10.1063/1.1515483 Optimization review
# [2] https://doi.org/10.1063/1.450914 Trust region method
# [3] 10.1007/978-0-387-40065-5 Numerical optimization
# [4] 10.1007/s00214-016-1847-3 Explorations of some refinements


import numpy as np

from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer


class RFOptimizer(HessianOptimizer):

    def optimize(self):
        gradient = self.geometry.gradient
        self.forces.append(-self.geometry.gradient)
        self.energies.append(self.geometry.energy)

        if self.cur_cycle > 0:
            self.update_trust_radius()
            self.update_hessian()

        H = self.H
        if self.geometry.internal:
            H_proj = self.geometry.internal.project_hessian(self.H)
            # Symmetrize hessian, as the projection probably breaks it?!
            H = (H_proj + H_proj.T) / 2
        # TODO: Neglect small eigenvalues with cartesian coordinates

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
        self.log(f"norm(step)={np.linalg.norm(step):.6f}")

        predicted_change = step.dot(gradient) + 0.5 * step.dot(self.H).dot(step)
        self.predicted_energy_changes.append(predicted_change)

        return step
