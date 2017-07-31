#!/usr/bin/env python3

import numpy as np

from cos.ChainOfStates import ChainOfStates

class Optimizer:

    def __init__(self, geometry, fix_ends=False, **kwargs):
        self.geometry = geometry
        self.fix_ends = fix_ends

        # Setting some default values
        self.max_cycles = 50
        self.max_force_thresh = 0.05
        self.rms_force_thresh = 0.01

        self.max_step = 0.04
        self.rel_step_thresh = 1e-3

        assert(self.max_step > self.rel_step_thresh)

        self.is_cos = issubclass(type(self.geometry), ChainOfStates)
        self.is_zts = getattr(self.geometry, "reparametrize", None)

        # Overwrite default values if they are supplied as kwargs
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.cur_cycle = 0
        self.coords = list()

        self.energies = list()
        self.forces = list()
        self.steps = list()
        self.max_forces = list()
        self.rms_forces = list()
        self.step_changes = [0, ]

    def print_convergence(self):
        print("cycle: {:04d} max(force): {:.5f} rms(force): {:.5f} "
                "d(step): {:.5f}".format(
            self.cur_cycle, self.max_forces[-1], self.rms_forces[-1],
            self.step_changes[-1])
        )

    def check_convergence(self):
        # Only use forces perpendicular to the mep
        if self.is_cos:
            forces = self.geometry.perpendicular_forces
        else:
            forces = self.forces[-1]

        max_force = forces.max()
        rms_force = np.sqrt(np.mean(np.square(forces)))

        self.max_forces.append(max_force)
        self.rms_forces.append(rms_force)

        self.is_converged = ((max_force <= self.max_force_thresh) and
                             (rms_force <= self.rms_force_thresh)
        )

    def check_step_change(self):
        if self.cur_cycle == 0:
            return

        step_change = np.linalg.norm(self.steps[-1] - self.steps[-2])
        self.step_changes.append(step_change)
        self.is_converged = step_change < self.rel_step_thresh

    def scale_by_max_step(self, steps):
        steps_max = steps.max()
        if steps_max > self.max_step:
            steps *= self.max_step / steps_max
        return steps

    def optimize(self):
        raise Exception("Not implemented!")

    def run(self):
        while True:
            if self.cur_cycle == self.max_cycles:
                print("Number of cycles exceeded!")
                break
            self.coords.append(self.geometry.coords)
            self.energies.append(self.geometry.energy)
            self.forces.append(self.geometry.forces)

            self.check_convergence()
            if self.is_converged:
                print("Converged! Gradients below threshold.")
                self.print_convergence()
                break

            steps = self.optimize()
            self.steps.append(steps)
            new_coords = self.geometry.coords + steps
            self.geometry.coords = new_coords

            if self.is_zts:
                self.geometry.reparametrize()

            self.check_step_change()
            if self.is_converged:
                print("Converged! Rel. step change below threshold.")
                self.print_convergence()
                break

            self.print_convergence()

            self.cur_cycle += 1
