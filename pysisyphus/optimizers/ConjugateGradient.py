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
        last_forces = self.forces[-1]
        if self.cur_cycle > 0:
            beta = (last_forces.dot(last_forces) /
                    self.forces[-2].dot(self.forces[-2])
            )
            steps = last_forces + beta*self.steps[-1]
            if np.isinf(beta):
                beta = 1.0
            steps = last_forces + beta*self.steps[-1]
        else:
            steps = last_forces
        steps = self.alpha * steps
        steps = self.scale_by_max_step(steps)

        last_coords = self.coords[-1]
        new_coords = last_coords + steps
        self.geometry.coords = new_coords

        new_forces = self.geometry.forces
        self.forces.append(new_forces)
        skip = self.backtrack(new_forces, last_forces)
        if skip:
            return None

        if self.align and self.is_cos:
            self.forces[-1], self.forces[-2], steps = self.fit_rigid((new_forces, last_forces, steps))
        #print(steps)
        return steps
