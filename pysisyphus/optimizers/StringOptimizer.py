#!/usr/bin/env python3

# [1] https://aip.scitation.org/doi/abs/10.1063/1.3664901
#     Behn, 2011, Freezing string method
# [2] https://aip.scitation.org/doi/pdf/10.1063/1.4804162
#     Zimmerman, 2013, Growing string with interpolation and optimization
#                      in internal coordiantes

import numpy as np
import scipy as sp

from pysisyphus.helpers import procrustes
from pysisyphus.optimizers.Optimizer import Optimizer

class StringOptimizer(Optimizer):

    def __init__(self, geometry, gamma=1.25, max_step=0.1, **kwargs):
        super().__init__(geometry, max_step=max_step, **kwargs)

        # gamma = 1.25 Hartree/Bohr² ~ 5 Hartree/Angstrom²
        self.gamma = gamma

    def prepare_opt(self):
        if self.is_cos and self.align:
            procrustes(self.geometry)

    def reset(self):
        pass

    def fit_rigid(self, vectors=None):
        rot_mats = procrustes(self.geometry)
        G = sp.linalg.block_diag(*rot_mats)
        rotated_vectors = [vec.dot(G) for vec in vectors]

        return rotated_vectors

    def restrict_step_components(self, steps):
        too_big = np.abs(steps) > self.max_step
        self.log(f"Found {np.sum(too_big)} big step components.")
        signs = np.sign(steps[too_big])
        steps[too_big] = signs * self.max_step
        return steps

    def optimize(self):
        forces = self.geometry.forces
        self.energies.append(self.geometry.energy)
        self.forces.append(forces)

        # Raises IndexError in cycle 0 and evaluates to False when the
        # string size changed after the last cycle.
        try:
            string_size_changed = forces.size != self.forces[-2].size
        except IndexError:
            string_size_changed = True
        # string_size_changed = True


        # if not string_size_changed and self.is_cos and self.align:
            # prev_forces = self.forces[-1]
            # prev_step = self.steps[-1]
            # vectors = (forces, prev_forces, prev_step)
            # forces, prev_forces, prev_step = self.fit_rigid(vectors)
            # self.forces[-1] = prev_forces
            # self.steps[-1] = prev_step
            # self.log("Aligned geometries")
        # self.forces.append(forces)

        sd_step = forces / self.gamma
        # Steepest descent in the first cbeginning
        # if (self.cur_cycle == 0) or string_size_changed:
        if True:
            step = sd_step
            self.log("Taking steepest descent step.")
        # Conjugate Gradient later one
        else:
            cur_norm = np.linalg.norm(forces)
            prev_norm = np.linalg.norm(self.forces[-2])
            quot = min(cur_norm**2 / prev_norm**2, 1)
            step = sd_step + quot*self.steps[-1]
            self.log("Taking conjugate gradient step.")

        step = self.restrict_step_components(step)

        return step
