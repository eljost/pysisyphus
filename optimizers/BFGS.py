#!/usr/bin/env python3

import numpy as np

from pysisyphus.optimizers.BacktrackingOptimizer import BacktrackingOptimizer

class BFGS(BacktrackingOptimizer):

    def __init__(self, geometry, **kwargs):
        super(BFGS, self).__init__(geometry, **kwargs)

        self.inv_hessian = np.eye(self.geometry.coords.size)

    def optimize(self):
        # Calculate initial forces in the first iteration
        if len(self.forces) is 0:
            last_forces = self.geometry.forces
        else:
            last_forces = self.forces[-1]

        steps = np.dot(self.inv_hessian, last_forces)
        steps = self.scale_by_max_step(steps)

        new_coords = self.coords + self.alpha * steps

        #coords_hold = self.coords.copy()
        #forces_hold = forces[-1].copy()

        self.forces.append(self.geometry.forces)
        self.energies.append(self.geometry.energy)
        forces = self.forces[-1]


        # !!!!!!!!!
        # fit_rigid
        # !!!!!!!!!

        #import sys; sys.exit()


        return steps
