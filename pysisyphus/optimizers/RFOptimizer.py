#!/usr/bin/env python3

# [1] http://aip.scitation.org/doi/10.1063/1.1515483

import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer

class RFOptimizer(Optimizer):

    def __init__(self, geometry, **kwargs):
        super().__init__(geometry, **kwargs)

    def prepare_opt(self):
        self.H = self.geometry.get_initial_hessian()

    def bfgs_update(self):
        # Eq. (44) in [1]
        dx = self.coords[-1] - self.coords[-2]
        dg = -(self.forces[-1] - self.forces[-2])
        second_term = np.outer(dg, dg) / np.inner(dg, dx)
        third_term = (self.H.dot(np.outer(dx, dx)).dot(self.H)
                      / dx.dot(self.H.dot(dx)))
        self.H += second_term - third_term

    def optimize(self):
        gradient = self.geometry.gradient
        self.forces.append(-self.geometry.gradient)
        self.energies.append(self.geometry.energy)

        # Update hessian
        if self.cur_cycle > 0:
            self.bfgs_update()

        # Eq. (56) in [1]
        aug_hess = np.bmat(
                    ((self.H, gradient[:,None]),
                     (gradient[None,:], [[0]]))
        )
        np.testing.assert_allclose(aug_hess, aug_hess.T)
        eigvals, eigvecs = np.linalg.eigh(aug_hess)
        # Select eigenvector corresponding to smallest eigenvalue
        # As the eigenvalues are sorted in ascending order eigvals.argmin()
        # should always give 0...
        aug_step = eigvecs[:,eigvals.argmin()]
        # aug_step is currently a matrix. Convert it to an array.
        aug_step = np.asarray(aug_step).flatten()
        # Scale aug_step so the last element equals 1
        aug_step /= aug_step[-1]
        step = self.scale_by_max_step(aug_step[:-1])
        #step = aug_step[:-1]
        return step
