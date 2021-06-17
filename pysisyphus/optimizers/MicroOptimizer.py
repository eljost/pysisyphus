import logging
import os
from pathlib import Path

import numpy as np

from pysisyphus.optimizers.closures import bfgs_multiply
from pysisyphus.optimizers.hessian_updates import double_damp
from pysisyphus.optimizers.poly_fit import poly_line_search


class MicroOptimizer:
    def __init__(
        self,
        geom,
        step="lbfgs",
        line_search=True,
        max_cycles=100_000_000,
        max_step=0.2,
        keep_last=10,
        rms_force=None,
        double_damp=True,
        dump=False,
        **kwargs,
    ):
        self.geometry = geom
        self.step_funcs = {
            "sd": self.sd_step,
            "cg": self.cg_step,
            "lbfgs": self.lbfgs_step,
        }
        self.step_func = self.step_funcs[step]
        self.line_search = line_search
        self.max_cycles = max_cycles
        self.max_step = max_step
        self.keep_last = keep_last
        self.double_damp = double_damp
        self.rms_force = rms_force
        self.dump = dump

        self.prev_energy = None
        self.prev_step = None
        self.prev_forces = None
        self.coord_diffs = list()
        self.grad_diffs = list()

        self.trj_fn = "microopt.trj"
        if Path(self.trj_fn).exists():
            os.remove(self.trj_fn)

        self.logger = logging.getLogger("optimizer")

        if kwargs:
            msg = "Got unsupported keyword arguments: " + ", ".join(
                [f"{k}={v}" for k, v in kwargs.items()]
            )
            self.log(msg)

    def log(self, msg):
        self.logger.debug(msg)

    def run(self):
        self.geometry.reparametrize()

        for self.cur_cycle in range(self.max_cycles):
            if self.dump:
                with open(self.trj_fn, "a") as handle:
                    handle.write(self.geometry.as_xyz() + "\n")

            results = self.geometry.get_energy_and_forces_at(self.geometry.coords)
            forces = results["forces"]
            energy = results["energy"]
            rms = np.sqrt(np.mean(forces ** 2))
            print(f"{self.cur_cycle:03d} rms(f)={rms:.6f}")
            if self.rms_force and rms <= self.rms_force:
                print("Converged!")
                self.is_converged = True
                break

            self.take_step(energy, forces)

    def take_step(self, energy, forces):
        self.log(
            f"Cycle {self.cur_cycle:03d}, energy={energy:.6f} au, "
            f"norm(forces)={np.linalg.norm(forces):.6f}"
        )

        ip_gradient, ip_step = None, None
        if self.line_search and (self.cur_cycle > 0):
            ip_energy, ip_gradient, ip_step = poly_line_search(
                energy,
                self.prev_energy,
                -forces,
                -self.prev_forces,
                self.prev_step,
            )
        # Use the interpolated gradient for the step if interpolation succeeded
        if (ip_gradient is not None) and (ip_step is not None):
            forces = -ip_gradient
        # Keep the original gradient when the interpolation failed, but use
        # zero (interpolated) step.
        else:
            ip_step = np.zeros_like(forces)

        # Calculate actual step with potentially interpolated forces
        if self.cur_cycle == 0:
            step = self.sd_step(forces)
        else:
            step = self.step_func(forces)

        # Form full step. If we did not interpolate or interpolation failed ip_step will be zero.
        step = step + ip_step
        # Restrict absolute value of step vector entries
        step *= min(self.max_step / np.abs(step).max(), 1)

        self.geometry.reparametrize()

        self.prev_step = step
        self.prev_energy = energy
        self.prev_forces = forces
        new_coords = self.geometry.coords + step
        self.geometry.coords = new_coords

    def sd_step(self, forces):
        step = forces
        return step

    def cg_step(self, forces):
        beta = forces.dot(forces) / self.prev_forces.dot(self.prev_forces)
        step = forces + beta * self.prev_step
        return step

    def lbfgs_step(self, forces):
        y = self.prev_forces - forces
        s = self.prev_step

        if self.double_damp:
            s, y = double_damp(s, y, s_list=self.coord_diffs, y_list=self.grad_diffs)
        self.grad_diffs.append(y)
        self.coord_diffs.append(s)

        # Drop superfluous oldest vectors
        self.coord_diffs = self.coord_diffs[-self.keep_last :]
        self.grad_diffs = self.grad_diffs[-self.keep_last :]

        step = bfgs_multiply(
            self.coord_diffs,
            self.grad_diffs,
            forces,
            gamma_mult=False,
            logger=self.logger,
        )

        return step
