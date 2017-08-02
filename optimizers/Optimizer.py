#!/usr/bin/env python3

import time

import numpy as np
import scipy as sp

from pysisyphus.cos.ChainOfStates import ChainOfStates
from qchelper.geometry import make_trj_str


gauss_loose = {
    "max_force_thresh": 2.5e-3,
    "rms_force_thresh": 1.7e-3,
    "max_step_thresh": 1.0e-2,
    "rms_step_thresh": 6.7e-3
}

class Optimizer:

    def __init__(self, geometry, convergence=gauss_loose, **kwargs):
        self.geometry = geometry

        self.convergence = convergence
        for key, value in convergence.items():
            setattr(self, key, value)

        # Setting some default values
        self.max_cycles = 50
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
        self.cycle_times = list()

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

        keys = self.convergence.keys()
        this_cycle = {
            "max_force_thresh": max_force,
            "rms_force_thresh": rms_force,
            "max_step_thresh": max_step,
            "rms_step_thresh": rms_step
        }

        self.is_converged = all(
            [this_cycle[key] <= getattr(self, key)
             for key in self.convergence.keys()
            ]
        )

    def print_header(self):
        hs = "max(force) rms(force) max(step) rms(step) s/cycle".split()
        header = "cycle" + " ".join([h.rjust(13) for h in hs])
        print(header)

    def print_convergence(self):
        int_fmt = "{:>5d}"
        float_fmt = "{:>12.6f}"
        conv_str = int_fmt + " " + (float_fmt + " ") * 4 + "{:>12.1f}"
        print(conv_str.format(
            self.cur_cycle, self.max_forces[-1], self.rms_forces[-1],
            self.max_steps[-1], self.rms_steps[-1], self.cycle_times[-1])
        )

    def scale_by_max_step(self, steps):
        steps_max = steps.max()
        if steps_max > self.max_step:
            steps *= self.max_step / steps_max
        return steps

    def procrustes(self):
        # https://github.com/pycogent/pycogent/blob/master/cogent/cluster/procrustes.py
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.procrustes.html
        # Center of geometry
        first_coords = self.geometry.images[0].coords
        first_centroid = np.mean(first_coords)
        first_coords_centered = first_coords - first_centroid
        print(first_coords_centered)

        for i in range(0, len(self.geometry.images)-1):
            ith_coords = self.geometry.images[i].coords
            ith_centroid = np.mean(ith_coords)
            ith_coords_centered = ith_coords - ith_centroid
            matrix = np.dot(ith_coords_centered.transpose(), first_coords_centered)

            U, W, Vt = sp.linalg.svd(matrix)
            print(i)

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
        self.print_header()
        while True:
            start_time = time.time()
            if self.cur_cycle == self.max_cycles:
                print("Number of cycles exceeded!")
                break

            self.coords.append(self.geometry.coords)

            steps = self.optimize()
            self.steps.append(steps)
            self.energies.append(self.geometry.energy)

            self.check_convergence()

            new_coords = self.geometry.coords + steps
            self.geometry.coords = new_coords

            self.cur_cycle += 1
            end_time = time.time()
            elapsed_seconds = end_time - start_time
            self.cycle_times.append(elapsed_seconds)

            self.save_cycle()
            self.print_convergence()
            if self.is_converged:
                print("Converged!")
                print()
                break

            if self.is_zts:
                self.geometry.reparametrize()

