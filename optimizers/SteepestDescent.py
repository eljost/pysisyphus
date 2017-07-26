#!/usr/bin/env python3

import numpy as np

from optimizers.Optimizer import Optimizer

class SteepestDescent(Optimizer):

    def __init__(self, geometry, **kwargs):
        self.alpha = -0.05
        super(SteepestDescent, self).__init__(geometry, **kwargs)

        assert(self.alpha < 0), "Alpha should be negative!"

    def optimize(self):
        if self.cur_cycle > 0:
            self.skip = self.backtrack()
        return self.alpha*self.forces[-1]
