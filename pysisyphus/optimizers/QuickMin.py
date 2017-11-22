#!/usr/bin/env python3

import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer

class QuickMin(Optimizer):

    def __init__(self, geometry, **kwargs):
        super(QuickMin, self).__init__(geometry, **kwargs)

        self.dt = 0.001

    def prepare_opt(self):
        self.velocities = [np.zeros_like(self.geometry.coords), ]

    def optimize(self):
        cur_velocities = self.velocities[-1]
        if self.align and self.is_cos:
            cur_velocities = self.fit_rigid((cur_velocities, ))
        cur_forces = self.geometry.forces
        self.forces.append(cur_forces)
        if self.cur_cycle == 0:
            tmp_velocities = np.zeros_like(cur_velocities)
        else:
            last_velocities = self.velocities[-2]
            coeff = last_velocities.dot(cur_forces)
            if coeff > 0:
                tmp_velocities = coeff*cur_forces/np.linalg.norm(cur_forces)
            else:
                tmp_velocities = np.zeros_like(cur_velocities)
        accelerations = cur_forces / self.geometry.masses_rep
        new_velocities = tmp_velocities + self.dt*accelerations
        steps = new_velocities*self.dt + 1/2*accelerations*self.dt**2
        steps = self.scale_by_max_step(steps)
        new_coords = self.geometry.coords + steps
        self.geometry.coords = new_coords
        self.velocities.append(new_velocities)

        return steps
