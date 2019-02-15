#!/usr/bin/env python3

import numpy as np

from pysisyphus.helpers import fit_rigid, procrustes
from pysisyphus.optimizers.BacktrackingOptimizer import BacktrackingOptimizer

# http://ikuz.eu/2015/04/15/the-concept-of-conjugate-gradient-descent-in-python/

class ConjugateGradient(BacktrackingOptimizer):

    def __init__(self, geometry, alpha=0.1, formula="FR", dont_skip=True,
                 **kwargs):
        super(ConjugateGradient, self).__init__(geometry, alpha=alpha,
                                                **kwargs)

        self.formula = formula
        self.dont_skip = dont_skip

    def prepare_opt(self):
        if self.is_cos and self.align:
            procrustes(self.geometry)
        # Calculate initial forces before the first iteration
        self.coords.append(self.geometry.coords)
        self.forces.append(self.geometry.forces)

    def reset(self):
        super().reset()
        # Check if the number of images changed
        if self.forces[-1].shape != self.coords[-1].shape:
            new_forces = self.geometry.forces
            self.forces.append(new_forces)
        self.resetted = True

    def get_beta(self, cur_forces, prev_forces):
        # Fletcher-Reeves formula
        if self.formula == "FR":
            beta = cur_forces.dot(cur_forces) / prev_forces.dot(prev_forces)
        # Polak-Ribiere
        elif self.formula == "PR":
            beta = (-cur_forces.dot(prev_forces-cur_forces)
                    / prev_forces.dot(prev_forces))
            beta_old = (cur_forces.dot(cur_forces-prev_forces)
                    / prev_forces.dot(prev_forces))
            self.log(f"beta_old={beta_old:.4f}, beta={beta:.4f}")
            if beta < 0:
                self.log(f"beta = {beta:.04f} < 0, resetting to 0")
                # beta = 0 basically restarts CG, as no previous step
                # information is mixed into the current step.
                beta = 0
        return beta

    def optimize(self):
        cur_forces = self.forces[-1]

        if not self.resetted and self.cur_cycle > 0:
            prev_forces = self.forces[-2]
            beta = self.get_beta(cur_forces, prev_forces)
            self.log(f"beta = {beta:.06f}")
            if np.isinf(beta):
                beta = 1.0
            steps = cur_forces + beta*self.steps[-1]
        else:
            # Start with steepest descent in the first iteration
            steps = cur_forces
            self.resetted = False
        steps = self.alpha * steps
        steps = self.scale_by_max_step(steps)

        last_coords = self.coords[-1]
        new_coords = last_coords + steps
        self.geometry.coords = new_coords

        new_forces = self.geometry.forces
        new_energy = self.geometry.energy

        skip = self.backtrack(new_forces, cur_forces)
        # Imho backtracking gives bad results here, so only use it if
        # explicitly requested (self.dont_skip == False).
        if (not self.dont_skip) and skip:
            self.geometry.coords = last_coords
            return None

        if self.align and self.is_cos:
            (new_forces, cur_forces, steps), _, _ = fit_rigid(self.geometry,
                                                              (new_forces,
                                                               cur_forces,
                                                               steps))
            self.geometry.coords -= steps
            # Set the calculated properties on the rotated geometries
            self.geometry.energy = new_energy
            self.geometry.forces = new_forces

        self.forces[-1] = cur_forces
        self.forces.append(new_forces)
        self.energies.append(new_energy)

        return steps
