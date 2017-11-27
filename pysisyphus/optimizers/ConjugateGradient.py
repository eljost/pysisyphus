#!/usr/bin/env python3

import numpy as np

from pysisyphus.optimizers.BacktrackingOptimizer import BacktrackingOptimizer

# http://ikuz.eu/2015/04/15/the-concept-of-conjugate-gradient-descent-in-python/

class ConjugateGradient(BacktrackingOptimizer):

    def __init__(self, geometry, **kwargs):
        super(ConjugateGradient, self).__init__(geometry, alpha=0.1, **kwargs)

    def prepare_opt(self):
        if self.is_cos and self.align:
            self.procrustes()
        # Calculate initial forces before the first iteration
        self.coords.append(self.geometry.coords)
        self.forces.append(self.geometry.forces)

    def optimize(self):
        cur_forces = self.forces[-1]
        if self.cur_cycle > 0:
            prev_forces = self.forces[-2]
            beta = cur_forces.dot(cur_forces) / prev_forces.dot(prev_forces)
            if np.isinf(beta):
                beta = 1.0
            steps = cur_forces + beta*self.steps[-1]
        else:
            steps = cur_forces
        steps = self.alpha * steps
        steps = self.scale_by_max_step(steps)

        last_coords = self.coords[-1]
        new_coords = last_coords + steps
        self.geometry.coords = new_coords

        new_forces = self.geometry.forces
        new_energy = self.geometry.energy

        skip = self.backtrack(new_forces, cur_forces)
        # Imho backtracking gives bad results here, so we don't use it.
        #if skip:
        #    return None

        if self.align and self.is_cos:
            new_forces, cur_forces, steps = self.fit_rigid((new_forces,
                                                            cur_forces, steps))
            self.geometry.coords -= steps
            # Set the calculated properties on the rotated geometries
            self.geometry.energy = new_energy
            self.geometry.forces = new_forces

        self.forces[-1] = cur_forces
        self.forces.append(new_forces)
        self.energies.append(new_energy)

        return steps
