#!/usr/bin/env python3

import numpy as np

from pysisyphus.optimizers.BacktrackingOptimizer import BacktrackingOptimizer

class SteepestDescent(BacktrackingOptimizer):

    def __init__(self, geometry, **kwargs):
        super(SteepestDescent, self).__init__(geometry, **kwargs)

    def optimize(self):
        # !!!!!!!!!
        # procrustes
        # !!!!!!!!!

        self.forces.append(self.geometry.forces)

        if self.cur_cycle > 0:
            self.skip = self.backtrack(self.forces[-1], self.forces[-2])

        step = self.alpha*self.forces[-1]
        step = self.scale_by_max_step(step)
        return step
