import logging
import os
from pathlib import Path

import numpy as np

from pysisyphus.optimizers.poly_fit import poly_line_search


class MicroOptimizer:
    def __init__(
        self,
        geom,
        step="sd",
        line_search=True,
        max_cycles=100_000_000,
        max_step=0.2,
        rms_force=None,
        dump=False,
        **kwargs,
    ):
        self.geometry = geom
        self.step_funcs = {
            "sd": self.sd_step,
            "cg": self.cg_step,
        }
        self.step_func = self.step_funcs[step]
        self.line_search = line_search
        self.max_cycles = max_cycles
        self.max_step = max_step
        self.rms_force = rms_force
        self.dump = dump

        self.prev_step = None
        self.prev_forces = None

        self.trj_fn = "ipiopt.trj"
        if Path(self.trj_fn).exists():
            os.remove(self.trj_fn)

        self.logger = logging.getLogger("optimizer")

    def log(self, msg):
        self.logger.debug(msg)

    def run(self):
        for self.cur_cycle in range(self.max_cycles):
            if self.dump:
                with open(self.trj_fn, "a") as handle:
                    handle.write(self.geometry.as_xyz() + "\n")

            results = self.geometry.get_energy_and_forces_at(self.geometry.coords)
            forces = results["forces"]
            energy = results["energy"]
            if self.rms_force and np.sqrt(np.mean(forces ** 2)) <= self.rms_force:
                print("Converged!")
                break
            rms = np.sqrt(np.mean(forces ** 2))
            print(f"{self.cur_cycle:03d} rms(f)={rms:.6f}")

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
        if (self.prev_step is None) or (self.prev_forces is None):
            step = self.sd_step(forces)
        else:
            beta = forces.dot(forces) / self.prev_forces.dot(self.prev_forces)
            step = forces + beta * self.prev_step
        return step
