#!/usr/bin/env python3

import logging

import os
from pathlib import Path
import time

import numpy as np
import pandas as pd
import scipy.linalg
import yaml

from pysisyphus.cos.ChainOfStates import ChainOfStates


gauss_loose = {
    "max_force_thresh": 2.5e-3,
    "rms_force_thresh": 1.7e-3,
    "max_step_thresh": 1.0e-2,
    "rms_step_thresh": 6.7e-3
}


class Optimizer:

    def __init__(self, geometry, convergence=gauss_loose,
                 align=False, dump=False,
                 climb=False, climb_multiple=5.0, climb_rms=False,
                 last_cycle=None,
                 **kwargs):
        self.geometry = geometry

        self.convergence = convergence
        self.align = align
        self.dump = dump
        self.climb = climb
        self.climb_multiple = climb_multiple
        self.climb_rms = climb_rms
        self.last_cycle = last_cycle

        for key, value in convergence.items():
            setattr(self, key, value)

        # Setting some default values
        self.started_climbing = False
        self.max_cycles = 50
        self.max_step = 0.04
        self.rel_step_thresh = 1e-3
        self.out_dir = os.getcwd()

        assert(self.max_step > self.rel_step_thresh)

        self.is_cos = issubclass(type(self.geometry), ChainOfStates)
        self.is_zts = getattr(self.geometry, "reparametrize", None)

        self.image_num = 1
        if self.is_cos:
            self.image_num = len(self.geometry.moving_indices)
            print(f"Path with {self.image_num} moving images.")

        # Overwrite default values if they are supplied as kwargs
        for key, value in kwargs.items():
            setattr(self, key, value)

        self.out_dir = Path(self.out_dir)
        if not self.out_dir.exists():
            os.mkdir(self.out_dir)
        self.logger = logging.getLogger("optimizer")

        self.list_attrs = "coords energies forces steps " \
                          "max_forces rms_forces max_steps rms_steps " \
                          "cycle_times tangents".split()
        for la in self.list_attrs:
            setattr(self, la, list())

        self.opt_results_fn = "optimizer_results.yaml"
        self.image_results_fn = "image_results.yaml"
        self.image_results = list()

        self.cur_cycle = 0
        if self.last_cycle:
            # Increase cycle number by one as we don't want to
            # redo the last cycle.
            self.cur_cycle = last_cycle + 1
            self.restart()

    def save_also(self):
        pass

    def restart(self):
        with open(self.image_results_fn) as handle:
            self.image_results = yaml.load(handle.read())
        with open(self.opt_results_fn) as handle:
            yaml_str = handle.read()
        opt_results = yaml.load(yaml_str)
        for key, val in opt_results.items():
            setattr(self, key, val)
        self.check_for_climbing_start()

    def log(self, message):
        self.logger.debug(f"Cycle {self.cur_cycle:03d}, {message}")

    def check_convergence(self, multiple=1.0):
        """Check if the current convergence of the optimization
        is equal to or below the required thresholds, or a multiple
        thereof. The latter is used in initiating the climbing image.
        """
        # When using a ChainOfStates method we are only interested
        # in optimizing the forces perpendicular to the MEP.
        if self.is_cos:
            forces = self.geometry.perpendicular_forces
        else:
            forces = self.forces[-1]
        step = self.steps[-1]

        max_force = np.abs(forces).max()
        rms_force = np.sqrt(np.mean(np.square(forces)))
        self.max_forces.append(max_force)
        self.rms_forces.append(rms_force)

        max_step = step.max()
        rms_step = np.sqrt(np.mean(np.square(step)))
        self.max_steps.append(max_step)
        self.rms_steps.append(rms_step)

        this_cycle = {
            "max_force_thresh": max_force,
            "rms_force_thresh": rms_force,
            "max_step_thresh": max_step,
            "rms_step_thresh": rms_step
        }

        return all(
            [this_cycle[key] <= getattr(self, key)*multiple
             for key in self.convergence.keys()]
        )

    def check_for_climbing_start(self):
        # Return False if we don't want to climb or are already
        # climbing.
        if (not self.climb) or self.started_climbing:
            return False

        # Only initiate climbing on a sufficiently converged MEP,
        # determined by a supplied threshold for the rms_force,
        # or as a multiple of all given convergence thresholds.
        #
        # If climb_rms was given we will use it for the check. Otherwise
        # we use a multiple of all four convergence thresholds.
        if self.climb_rms:
            climb_now = self.rms_forces[-1] <= self.climb_rms
        else:
            climb_now = self.check_convergence(self.climb_multiple)

        self.started_climbing = climb_now
        self.geometry.climb = climb_now
        if climb_now:
            self.log("starting to climb in next iteration.")
            print("starting to climb in next iteration.")
        return climb_now

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
        steps_max = np.abs(steps).max()
        if steps_max > self.max_step:
            steps *= self.max_step / steps_max
        return steps

    def prepare_opt(self):
        pass

    def optimize(self):
        raise Exception("Not implemented!")

    def write_to_out_dir(self, out_fn, content, mode="w"):
        out_path = self.out_dir / out_fn
        with open(out_path, mode) as handle:
            handle.write(content)

    def write_image_trjs(self):
        base_name = "image_{:03d}.trj"
        for i, image in enumerate(self.geometry.images):
            image_fn = base_name.format(i)
            comment = f"cycle {self.cur_cycle}"
            as_xyz = image.as_xyz(comment)
            self.write_to_out_dir(image_fn, as_xyz+"\n", "a")

    def write_results(self):
        # Save results from the Geometry.
        self.image_results.append(self.geometry.results)
        self.write_to_out_dir(self.image_results_fn,
                              yaml.dump(self.image_results))

        # Save results from the Optimizer
        opt_results = {la: getattr(self, la) for la in self.list_attrs}
        opt_results.update(self.save_also())
        self.write_to_out_dir(self.opt_results_fn,
                              yaml.dump(opt_results))

    def write_cycle_to_file(self):
        as_xyz_str = self.geometry.as_xyz()

        if self.is_cos:
            out_fn = "cycle_{:03d}.trj".format(self.cur_cycle)
            self.write_to_out_dir(out_fn, as_xyz_str)
            # Also write separate .trj files for every image in the cos
            self.write_image_trjs()
            self.write_results()
        else:
            # Append to .trj file
            out_fn = "optimization.trj"
            self.write_to_out_dir(out_fn, as_xyz_str+"\n", mode="a")

    def run(self):
        if not self.last_cycle:
            prep_start_time = time.time()
            self.prepare_opt()
            prep_end_time = time.time()
            prep_time = prep_end_time - prep_start_time
            print(f"Spent {prep_time:.1f} s preparing the first cycle.")

        self.print_header()
        while True:
            start_time = time.time()
            if self.cur_cycle == self.max_cycles:
                print("Number of cycles exceeded!")
                break

            self.coords.append(self.geometry.coords)

            steps = self.optimize()

            if steps is None:
                # Remove the previously added coords
                self.coords.pop(-1)
                continue

            if self.is_cos:
                self.tangents.append(self.geometry.get_tangents())

            self.steps.append(steps)

            # Convergence check
            self.is_converged = self.check_convergence()

            self.check_for_climbing_start()
            end_time = time.time()
            elapsed_seconds = end_time - start_time
            self.cycle_times.append(elapsed_seconds)

            if self.dump:
                self.write_cycle_to_file()

            self.print_convergence()
            if self.is_converged:
                print("Converged!")
                print()
                break

            new_coords = self.geometry.coords + steps
            self.geometry.coords = new_coords

            if self.is_zts:
                self.geometry.reparametrize()

            stop_signs = ("stop", "STOP")
            for ss in stop_signs:
                if os.path.exists(ss):
                    print("Found stop sign. Stopping optimization.")
                    os.remove(ss)
                    return

            self.cur_cycle += 1
