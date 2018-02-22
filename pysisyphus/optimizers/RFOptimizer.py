#!/usr/bin/env python3

# [1] http://aip.scitation.org/doi/10.1063/1.1515483 Optimization review
# [2] https://doi.org/10.1063/1.450914 Trust region method
# [3] 10.1007/978-0-387-40065-5 Numerical optimization

import logging

import numpy as np
from scipy.optimize import minimize

from pysisyphus.optimizers.Optimizer import Optimizer

class RFOptimizer(Optimizer):

    def __init__(self, geometry, trust_radius=0.3, **kwargs):
        super().__init__(geometry, **kwargs)

        self.trust_radius = trust_radius
        self.min_trust_radius = 0.25*self.trust_radius
        self.trust_radius_max = 5*self.trust_radius
        self.predicted_energies = list()

        self.hessians = list()
        self.trust_radii = list()
        self.rfo_steps = list()
        self.predicted_energies = list()#[self.quadratic_approx(self.geometry.coords), ]

    def prepare_opt(self):
        self.H = self.geometry.get_initial_hessian()

    def keep(self):
        self.hessians.append(self.H.copy())
        self.trust_radii.append(self.trust_radius)
        self.log("!Saving hessian every iteration!")

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
        step = res.x
        self.log(f"optimized step norm: {np.linalg.norm(step):.4f}")
        if not res.success:
            raise Exception("LQA optimization failed!")
        return step

    def update_trust_radius(self):
        # [3] Chapter 4, Algorithm 4.1
        actual_reduction = self.energies[-2] - self.energies[-1]
        predicted_reduction = (self.predicted_energies[-2]
                               - self.predicted_energies[-1])
        reduction_ratio = actual_reduction / predicted_reduction
        self.log(f"reduction_ratio: {reduction_ratio:.3f}")
        if reduction_ratio < 0.25:
            self.trust_radius = max(0.25*self.trust_radius,
                                    self.min_trust_radius)
        elif reduction_ratio > 0.5:
            self.trust_radius = min(2*self.trust_radius, self.trust_radius_max)
        #if reduction_ratio < 0:
        #    self.log("step rejected")
        self.log(f"trust_radius {self.trust_radius:.2f}")

    def optimize(self):
        gradient = self.geometry.gradient
        self.forces.append(-self.geometry.gradient)
        self.energies.append(self.geometry.energy)

        predicted_energy = self.quadratic_approx(self.coords[-1])
        self.predicted_energies.append(predicted_energy)

        if self.cur_cycle > 0:
            self.bfgs_update()
            self.update_trust_radius()

        if self.geometry.internal:
            self.H = self.geometry.internal.project_hessian(self.H)
            # Symmetrize hessian, as the projection probably breaks it.
            self.H = (self.H + self.H.T) / 2

        # Eq. (56) in [1]
        aug_hess = np.bmat(
                    ((self.H, gradient[:,None]),
                     (gradient[None,:], [[0]]))
        )
        #aug_hess = (aug_hess+aug_hess.T)/2
        np.testing.assert_allclose(aug_hess, aug_hess.T)
        eigvals, eigvecs = np.linalg.eigh(aug_hess)
        # Select eigenvector corresponding to smallest eigenvalue.
        # As the eigenvalues are sorted in ascending order eigvals.argmin()
        # should always give 0...
        aug_step = eigvecs[:,eigvals.argmin()]
        # aug_step is currently a matrix. Convert it to an array.
        aug_step = np.asarray(aug_step).flatten()
        # Scale aug_step so the last element equals 1
        aug_step /= aug_step[-1]
        step = aug_step[:-1]

        self.keep()
        self.rfo_steps.append(step)

        step_norm = np.linalg.norm(step)
        # We use a trust region method instead
        #step = self.scale_by_max_step(step)
        if step_norm > self.trust_radius:
            step = self.find_step(step)
        self.log(f"norm(step) {np.linalg.norm(step):.4f}")

        return step
