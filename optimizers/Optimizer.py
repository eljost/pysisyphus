#!/usr/bin/env python3

import numpy as np

from cos.ChainOfStates import ChainOfStates

class Optimizer:

    def __init__(self, geometry, **kwargs):
        self.geometry = geometry

        # Setting some default values
        self.max_cycles = 15
        self.max_force_thresh = 0.01
        self.rms_force_thresh = 0.001

        self.max_step = 0.04

        self.is_cos = issubclass(type(self.geometry), ChainOfStates)
        self.is_zts = getattr(self.geometry, "reparametrize", None)

        # Overwrite default values if they are supplied as kwargs
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.cur_cycle = 0
        self.coords = list()

        self.forces = list()
        self.steps = list()
        self.max_forces = list()
        self.rms_forces = list()
        self.step_changes = [0, ]

        # Check if geometry defines it's own convergence check, e.g.
        # as in CoS methods where we're interested in the perpendicular
        # component of the force along the MEP.
        geom_conv_check = getattr(self.geometry, "check_convergence", None)
        if geom_conv_check:
            self.check_convergence = geom_conv_check

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
        self.is_converged = step_change < 1e-4

    def scale_by_max_step(self, steps):
        steps_max = steps.max()
        if steps_max > self.max_step:
            steps *= self.max_step / steps_max
        return steps

    def optimize(self):
        raise Exception("Not implemented!")

    def run(self):
        while self.cur_cycle < self.max_cycles:
            forces = self.geometry.forces
            self.forces.append(forces)
            self.coords.append(self.geometry.coords)

            self.check_convergence()
            if self.is_converged:
                break

            steps = self.optimize()
            steps = self.scale_by_max_step(steps)
            self.steps.append(steps)
            new_coords = self.geometry.coords + steps
            self.geometry.coords = new_coords

            if self.is_zts:
                self.geometry.reparametrize()

            self.check_step_change()
            if self.is_converged:
                break

            self.print_convergence()

            self.cur_cycle += 1
