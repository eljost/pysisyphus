#!/usr/bin/env python3

# [1] http://aip.scitation.org/doi/10.1063/1.1515483 Optimization review
# [2] https://doi.org/10.1063/1.450914 Trust region method
# [3] 10.1007/978-0-387-40065-5 Numerical optimization
# [4] 10.1007/s00214-016-1847-3 Explorations of some refinements


import numpy as np
from scipy.optimize import minimize

from pysisyphus.optimizers.Optimizer import Optimizer


class RFOptimizer(Optimizer):

    def __init__(self, geometry, trust_radius=0.3, **kwargs):
        super().__init__(geometry, **kwargs)

        self.trust_radius = trust_radius
        self.min_trust_radius = 0.25*self.trust_radius
        self.trust_radius_max = 5*self.trust_radius
        self.predicted_energy_changes = list()
        self.actual_energy_changes = list()

        # self.hessians = list()
        self.trust_radii = list()
        self.rfo_steps = list()

    def prepare_opt(self):
        self.H = self.geometry.get_initial_hessian()

    def keep(self):
        # self.hessians.append(self.H.copy())
        self.trust_radii.append(self.trust_radius)
        # self.log("!Saving hessian every iteration!")

    def bfgs_update(self):
        # Eq. (44) in [1]
        dx = self.coords[-1] - self.coords[-2]
        dg = -(self.forces[-1] - self.forces[-2])
        second_term = np.outer(dg, dg) / np.inner(dg, dx)
        third_term = (self.H.dot(np.outer(dx, dx)).dot(self.H)
                      / dx.dot(self.H.dot(dx)))
        self.H += second_term - third_term

    def quadratic_approx(self, step):
        E0 = self.energies[-1]
        g = -self.forces[-1]
        return E0 + np.inner(g, step) + 0.5*step.dot(self.H.dot(step))

    def predicted_change(self, step):
        E0 = self.energies[-1]
        g = -self.forces[-1]
        return np.inner(g, step) + 0.5*step.dot(self.H.dot(step))

    def find_step(self, step_guess):
        ineq_fun = lambda step: self.trust_radius - np.linalg.norm(step)
        constr = {
                "type": "ineq",
                "fun": ineq_fun,
        }
        res = minimize(self.quadratic_approx,
                 step_guess,
                 method="SLSQP",
                 constraints=constr)
        if not res.success:
            self.log("LQA optimization failed!")
        step = res.x
        self.log(f"LQA minimum: {res.fun:.6f} au")
        self.log(f"Optimized step norm: {np.linalg.norm(step):.4f}")
        return step

    def update_trust_radius(self):
        # [3] Chapter 4, Algorithm 4.1
        actual_reduction = self.energies[-2] - self.energies[-1]
        predicted = self.predicted_energy_changes[-1]
        actual = self.actual_energy_changes[-1]
        reduction_ratio = actual / predicted
        self.log(f"Predicted energy reduction: {predicted:.4e} au")
        self.log(f"Actual energy reduction: {actual:.4e} au")
        self.log(f"Reduction ratio: {reduction_ratio:.5f}")
        last_step_norm = np.linalg.norm(self.steps[-1])
        if reduction_ratio < 0.25:
            self.trust_radius = max(0.25*self.trust_radius,
                                    self.min_trust_radius)
            self.log("Decreasing trust radius.")
        # Only increase trust radius if last step norm was at least 80% of it
        # See [4], Appendix, step size and direction control
        elif reduction_ratio > 0.5 and (last_step_norm >= .8*self.trust_radius):
            self.trust_radius = min(2*self.trust_radius, self.trust_radius_max)
            self.log("Increasing trust radius.")
        self.log(f"Trust radius: {self.trust_radius:.3f}")

    def optimize(self):
        gradient = self.geometry.gradient
        self.forces.append(-self.geometry.gradient)
        self.energies.append(self.geometry.energy)

        if self.cur_cycle > 0:
            # Predicted changes
            # predicted_energy_change = self.quadratic_approx(self.steps[-1])
            predicted_energy_change = self.predicted_change(self.steps[-1])
            self.predicted_energy_changes.append(predicted_energy_change)
            # Actual change
            self.log(f"Last two energies: {self.energies[-2:]}")
            actual_energy_change = self.energies[-1] - self.energies[-2]
            self.actual_energy_changes.append(actual_energy_change)
            self.update_trust_radius()
            self.bfgs_update()

        if self.geometry.internal:
            self.H = self.geometry.internal.project_hessian(self.H)
            # Symmetrize hessian, as the projection probably breaks it.
            self.H = (self.H + self.H.T) / 2

        # Eq. (56) in [1]
        aug_hess = np.bmat(
                    ((self.H, gradient[:,None]),
                     (gradient[None,:], [[0]]))
        )
        eigvals, eigvecs = np.linalg.eigh(aug_hess)
        self.log(f"First 5 eigenvalues: {np.array2string(eigvals[:5], precision=2, suppress_small=True)}")
        self.log(f"Number of negative eigenvalues: {eigvals[eigvals < 0].size}")
        self.log(f"Lowest eigenvalue: {eigvals[0]:.6f}")
        # Select eigenvector corresponding to smallest eigenvalue.
        # As the eigenvalues are sorted in ascending order eigvals.argmin()
        # should always give 0...
        assert eigvals.argmin() == 0
        aug_step = eigvecs[:,0]
        # aug_step is currently a matrix. Convert it to an array.
        aug_step = np.asarray(aug_step).flatten()
        # Scale aug_step so the last element equals 1
        lambda_ = aug_step[-1]
        self.log(f"lambda: {lambda_:.6f}")
        aug_step /= lambda_
        step = aug_step[:-1]

        self.keep()
        self.rfo_steps.append(step)

        # Restrict elements of the the step vector to an allowed maximum
        # if they exceed it.
        # step[np.abs(step) > 0.3] = 0.3
        step_norm = np.linalg.norm(step)
        self.log(f"Unscaled norm(step): {step_norm:.4f}")
        # We use a trust region method instead
        if step_norm > self.trust_radius:
            step = self.find_step(step)
        self.log(f"norm(step): {np.linalg.norm(step):.4f}")
        self.log("")

        return step
