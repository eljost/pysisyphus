#!/usr/bin/env python3

import numpy as np

class Optimizer:

    def __init__(self, geometry, les=15, **kwargs):
        self.geometry = geometry

        # Setting some default values
        self.max_cycles = 15
        self.max_force_thresh = 0.01
        self.rms_force_thresh = 0.001
        # Overwrite default values if they are supplied as kwargs
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.cur_cycle = 0
        self.coords = list()

        self.forces = list()
        self.steps = list()
        self.max_forces = list()
        self.rms_forces = list()

    def check_convergence(self, forces):
        max_force = forces.max()
        rms_force = np.sqrt(np.mean(np.square(forces)))

        self.max_forces.append(max_force)
        self.rms_forces.append(rms_force)

        print("cycle: {:04d} max(force): {:.5f} rms(force): {:.5f}".format(
            self.cur_cycle, max_force, rms_force)
        )

        return ((max_force <= self.max_force_thresh) and
                (rms_force <= self.rms_force_thresh)
        )

    def optimize(self):
        raise Exception("Not implemented!")

    def run(self):
        while self.cur_cycle < self.max_cycles:
            forces = self.geometry.forces
            self.forces.append(forces)
            self.coords.append(self.geometry.coords)
            if self.check_convergence(forces):
                break
            step = self.optimize()
            self.steps.append(step)
            new_coords = self.geometry.coords + step
            self.geometry.coords = new_coords

            self.cur_cycle += 1
