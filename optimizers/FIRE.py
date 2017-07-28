#!/usr/bin/env python3

import numpy as np

from optimizers.BacktrackingOptimizer import BacktrackingOptimizer

class FIRE(BacktrackingOptimizer):
    # https://doi.org/10.1103/PhysRevLett.97.170201

    def __init__(self, geometry, **kwargs):
        self.dt = 0.1
        self.dt_max = 1
        self.N_acc = 5
        self.f_inc = 1.1
        self.f_acc = 0.99
        self.f_dec = 0.5

        self.n_reset = 0
        self.a_start = 0.1
        self.a = self.a_start

        self.v = np.zeros_like(geometry.coords)
        self.velocities = [self.v, ]

        super(FIRE, self).__init__(geometry, **kwargs)

    def optimize(self):
        forces = self.forces[-1]
        mixed_velocity = (
            (1.0 - self.a) * self.v +
             self.a * np.sqrt(
                        np.vdot(self.v, self.v) / np.vdot(forces, forces)
                      ) * forces
        )
        last_velocity = self.velocities[-1]
        if (self.cur_cycle > 0) and (np.vdot(last_velocity, forces) > 0):
            if self.n_reset > self.N_acc:
                self.dt = min(self.dt * self.f_inc, self.dt_max)
                self.a *= self.f_acc
            self.n_reset += 1
        else:
            mixed_velocity = np.zeros_like(forces)
            self.a = self.a_start
            self.dt *= self.f_acc
            self.n_reset = 0
        velocity = mixed_velocity + self.dt * forces
        self.velocities.append(velocity)
        steps = self.dt * velocity
        steps = self.scale_by_max_step(steps)

        return steps
