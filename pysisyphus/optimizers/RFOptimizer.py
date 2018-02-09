#!/usr/bin/env python3

import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer

class RFOptimizer(Optimizer):

    def __init__(self, geometry, **kwargs):
        super().__init__(geometry, **kwargs)

    def prepare_opt(self):
        self.hessian = self.geometry.get_initial_hessian()

    #def scale_by_max_step(self, steps):
    #    return steps

    def optimize(self):
        gradient = self.geometry.gradient
        self.coords.append(self.geometry.coords)
        self.forces.append(-self.geometry.gradient)
        self.energies.append(self.geometry.energy)
        aug_hess = np.bmat(
                    ((self.hessian, gradient[:,None]),
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
        #import pdb; pdb.set_trace()
        return step

        pass
        #self.forces.append(self.geometry.forces)
        #self.energies.append(self.geometry.energy)

        #return step
