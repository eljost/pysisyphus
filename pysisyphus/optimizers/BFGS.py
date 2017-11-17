#!/usr/bin/env python3

import numpy as np

from pysisyphus.optimizers.BacktrackingOptimizer import BacktrackingOptimizer

class BFGS(BacktrackingOptimizer):

    def __init__(self, geometry, **kwargs):
        super(BFGS, self).__init__(geometry, alpha=1.0, **kwargs)

        self.reset_hessian()
        self.eye = self.inv_hessian.copy()

    def reset_hessian(self):
        self.inv_hessian = np.eye(self.geometry.coords.size)

    def prepare_opt(self):
        if self.is_cos and self.align:
            self.procrustes()
        # Calculate initial forces before the first iteration
        self.coords.append(self.geometry.coords)
        self.forces.append(self.geometry.forces)
        self.energies.append(self.geometry.energy)

    def optimize(self):
        last_coords = self.coords[-1]
        last_forces = self.forces[-1]
        last_energy = self.energies[-1]

        steps = self.inv_hessian.dot(last_forces)
        steps = self.scale_by_max_step(steps)
        steps *= self.alpha

        new_coords = last_coords + steps
        self.geometry.coords = new_coords

        if self.is_cos and self.align:
            (last_coords, last_forces), self.inv_hessian = self.fit_rigid((last_coords,
                                                                           last_forces),
                                                                           self.inv_hessian)

        new_forces = self.geometry.forces
        new_energy = self.geometry.energy
        skip = self.backtrack(new_forces, last_forces)

        if skip:
            self.reset_hessian()
            self.geometry.coords = last_coords
            return None
        else:
            self.geometry.coords = last_coords
            self.geometry.energy = new_energy

            self.forces.append(new_forces)
            sigma = new_coords - last_coords
            forces_diff = -new_forces - (-last_forces)
            rho = 1.0 / np.dot(forces_diff, sigma)
            if np.array_equal(self.inv_hessian, self.eye):
                self.inv_hessian = (np.dot(forces_diff, sigma) /
                                    np.dot(forces_diff, forces_diff) *
                                    self.eye
                )
            # Inverse hessian update
            A = (self.eye -
                 sigma[:,np.newaxis] * forces_diff[np.newaxis,:] * rho
            )
            B = (self.eye -
                 forces_diff[:,np.newaxis] * sigma[np.newaxis,:] * rho
            )
            self.inv_hessian = (
                    np.dot(A, np.dot(self.inv_hessian, B)) +
                    sigma[:,np.newaxis] * sigma[np.newaxis,:] * rho
            )

        return steps
