#!/usr/bin/env python3

import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer

class QuickMin(Optimizer):

    def __init__(self, geometry, dt=0.35, **kwargs):
        super(QuickMin, self).__init__(geometry, **kwargs)

        self.dt = dt

    def prepare_opt(self):
        self.velocities = [np.zeros_like(self.geometry.coords), ]

    def optimize(self):
        if self.align and self.is_cos:
            (self.velocities[-1], ) = self.fit_rigid((self.velocities[-1], ))

        cur_velocities = self.velocities[-1]
        cur_forces = self.geometry.forces
        self.forces.append(cur_forces)
        self.energies.append(self.geometry.energy)

        if self.cur_cycle == 0:
            tmp_velocities = np.zeros_like(cur_velocities)
        else:
            overlap = self.velocities[-2].dot(cur_forces)
            if overlap > 0:
                tmp_velocities = (overlap * cur_forces
                                  / cur_forces.dot(cur_forces))
            else:
                tmp_velocities = np.zeros_like(cur_velocities)
                self.log("resetted velocities")

        accelerations = cur_forces / self.geometry.masses_rep
        new_velocities = tmp_velocities + self.dt*accelerations
        steps = new_velocities*self.dt + 1/2*accelerations*self.dt**2
        steps = self.scale_by_max_step(steps)
        self.velocities.append(new_velocities)
        velo_norm = np.linalg.norm(new_velocities)
        acc_norm = np.linalg.norm(accelerations)
        self.log(f"norm(v) = {velo_norm:.4f}, norm(a) = {acc_norm:.4f}")

        return steps
