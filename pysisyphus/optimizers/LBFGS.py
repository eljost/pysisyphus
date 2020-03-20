#!/usr/bin/env python3

import numpy as np

from pysisyphus.helpers import fit_rigid
from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.closures import bfgs_multiply
from pysisyphus.optimizers.restrict_step import scale_by_max_step


class LBFGS(Optimizer):
    def __init__(self, geometry, alpha=1.0, keep_last=7, beta=1, max_step=0.2,
                  **kwargs):
        """[1] Nocedal, Wright - Numerical Optimization, 2006"""
        self.alpha = alpha
        self.beta = beta
        self.keep_last = int(keep_last)

        self.coord_diffs = list()
        self.grad_diffs = list()

        super().__init__(geometry, max_step=max_step, **kwargs)

    def _get_opt_restart_info(self):
        opt_restart_info = {
            "coord_diffs": np.array(self.coord_diffs).tolist(),
            "grad_diffs": np.array(self.grad_diffs).tolist()
        }
        return opt_restart_info

    def _set_opt_restart_info(self, opt_restart_info):
        self.coord_diffs = [np.array(cd) for cd in opt_restart_info["coord_diffs"]]
        self.grad_diffs = [np.array(gd) for gd in opt_restart_info["grad_diffs"]]

    def optimize(self):
        if self.is_cos and self.align:
            rot_vecs, rot_vec_lists, _ = fit_rigid(
                self.geometry,
                vector_lists=(self.steps, self.forces, self.coord_diffs, self.grad_diffs)
            )
            rot_steps, rot_forces, rot_coord_diffs, rot_grad_diffs = rot_vec_lists
            self.steps = rot_steps
            self.forces = rot_forces
            self.coord_diffs = rot_coord_diffs
            self.grad_diffs = rot_grad_diffs

        forces = self.geometry.forces
        self.forces.append(forces)
        energy = self.geometry.energy
        self.energies.append(energy)

        if self.cur_cycle > 0:
            prev_forces = self.forces[-2]
            grad_diff = -forces - -prev_forces
            self.grad_diffs.append(grad_diff)
            self.grad_diffs = self.grad_diffs[-self.keep_last:]

        step = -bfgs_multiply(self.coord_diffs, self.grad_diffs, forces, beta=self.beta)
        step = scale_by_max_step(step, self.max_step)

        self.coord_diffs.append(step)
        self.coord_diffs = self.coord_diffs[-self.keep_last:]

        return step
