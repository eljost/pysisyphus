#!/usr/bin/env python3

import numpy as np

from pysisyphus.helpers import fit_rigid, procrustes
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.closures import bfgs_multiply

# [1] Nocedal, Wright - Numerical Optimization, 2006

class LBFGS(Optimizer):
    def __init__(self, geometry, alpha=1.0, keep_last=15, **kwargs):
        self.alpha = alpha
        self.keep_last = keep_last
        super().__init__(geometry, **kwargs)

        self.sigmas = list()
        self.grad_diffs = list()

    def prepare_opt(self):
        if self.is_cos and self.align:
            procrustes(self.geometry)
        # Calculate initial forces before the first iteration
        self.coords.append(self.geometry.coords)
        self.forces.append(self.geometry.forces)
        self.energies.append(self.geometry.energy)

    def scale_by_max_step(self, steps):
        steps_max = np.abs(steps).max()
        if steps_max > self.max_step:
            fact = self.max_step / steps_max
            steps *= self.max_step / steps_max
        return steps

    def optimize(self):
        prev_coords = self.coords[-1]
        prev_forces = self.forces[-1]

        step = -bfgs_multiply(self.sigmas, self.grad_diffs, prev_forces)
        step = self.scale_by_max_step(step)

        new_coords = prev_coords + self.alpha*step

        coords_tmp = prev_coords.copy()
        forces_tmp = prev_forces.copy()

        self.geometry.coords = new_coords
        if self.is_cos and self.align:
            rot_vecs, rot_vec_lists, _ = fit_rigid(
                self.geometry,
                (prev_coords, prev_forces),
                vector_lists=(self.sigmas, self.grad_diffs)
            )
            prev_coords, prev_forces = rot_vecs
            self.sigmas, self.grad_diffs = rot_vec_lists

        new_forces = self.geometry.forces
        new_energy = self.geometry.energy

        sigma = new_coords - prev_coords
        self.sigmas.append(sigma)
        grad_diff = prev_forces - new_forces
        self.grad_diffs.append(grad_diff)

        self.sigmas = self.sigmas[-self.keep_last:]
        self.grad_diffs = self.grad_diffs[-self.keep_last:]

        # Because we add the step later on we restore the original
        # coordinates and set the appropriate energies and forces.
        self.geometry.coords = prev_coords
        self.geometry.forces = new_forces
        self.geometry.energy = new_energy

        self.forces.append(new_forces)
        self.energies.append(new_energy)

        return step

    def save_also(self):
        return {
            "alpha": self.alpha,
        }
