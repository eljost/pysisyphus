#!/usr/bin/env python3

import numpy as np

from cos.ChainOfStates import ChainOfStates
from qchelper.geometry import make_trj_str

class Optimizer:

    def __init__(self, geometry, fix_ends=False, **kwargs):
        self.geometry = geometry
        self.fix_ends = fix_ends

        # Setting some default values
        self.max_cycles = 50
        # Gaussian loose
        self.max_force_thresh = 2.5e-3
        self.rms_force_thresh = 1.7e-3
        self.max_step_thresh = 1.0e-2
        self.rms_step_thresh = 6.7e-3

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
        self.max_steps = list()
        self.rms_steps = list()

    def check_convergence(self):
        # Only use forces perpendicular to the mep
        if self.is_cos:
            forces = self.geometry.perpendicular_forces
        else:
            forces = self.forces[-1]
        step = self.steps[-1]

        max_force = forces.max()
        rms_force = np.sqrt(np.mean(np.square(forces)))
        self.max_forces.append(max_force)
        self.rms_forces.append(rms_force)

        max_step = step.max()
        rms_step = np.sqrt(np.mean(np.square(step)))
        self.max_steps.append(max_step)
        self.rms_steps.append(rms_step)

        self.is_converged = ((max_force <= self.max_force_thresh) and
                             (rms_force <= self.rms_force_thresh) and
                             (max_step <= self.max_step_thresh) and
                             (rms_step <= self.rms_step_thresh)
        )

    def print_convergence(self):
        print("cycle: {:04d} max(force): {:.5f} rms(force): {:.5f} "
                "max(step): {:.5f} rms(step): {:.5f}".format(
            self.cur_cycle, self.max_forces[-1], self.rms_forces[-1],
            self.max_steps[-1], self.rms_steps[-1])
        )

    def scale_by_max_step(self, steps):
        steps_max = steps.max()
        if steps_max > self.max_step:
            steps *= self.max_step / steps_max
        return steps

    def optimize(self):
        raise Exception("Not implemented!")

    def save_cycle(self):
        as_xyz_str = self.geometry.as_xyz()

        if self.is_cos:
            out_fn = "cycle_{:03d}.trj".format(self.cur_cycle)
            with open(out_fn, "w") as handle:
                handle.write(as_xyz_str)
        else:
            out_fn = "opt.trj"
            with open(out_fn, "a") as handle:
                handle.write(as_xyz_str)
                handle.write("\n")

    def run(self):
        while True:
            if self.cur_cycle == self.max_cycles:
                print("Number of cycles exceeded!")
                break

            self.coords.append(self.geometry.coords)
            self.forces.append(self.geometry.forces)
            self.energies.append(self.geometry.energy)


            steps = self.optimize()
            self.steps.append(steps)

            self.check_convergence()

            new_coords = self.geometry.coords + steps
            self.geometry.coords = new_coords

            self.save_cycle()
            self.print_convergence()
            if self.is_converged:
                print("Converged!")
                break

            if self.is_zts:
                self.geometry.reparametrize()

            self.cur_cycle += 1
