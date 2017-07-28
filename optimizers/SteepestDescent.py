#!/usr/bin/env python3

import numpy as np

from optimizers.BacktrackingOptimizer import BacktrackingOptimizer

class SteepestDescent(BacktrackingOptimizer):

    def __init__(self, geometry, **kwargs):
        self.alpha = 0.05
        super(SteepestDescent, self).__init__(geometry, **kwargs)

        assert(self.alpha > 0), "Alpha should be positive!"

    def optimize(self):
        if self.cur_cycle > 0:
            self.skip = self.backtrack()
        step = self.alpha*self.forces[-1]
        step = self.scale_by_max_step(step)
        return step
