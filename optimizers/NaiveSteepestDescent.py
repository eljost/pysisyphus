#!/usr/bin/env python3

from pysisyphus.optimizers.Optimizer import Optimizer

class NaiveSteepestDescent(Optimizer):

    def __init__(self, geometry, **kwargs):
        self.alpha = 0.05
        super(NaiveSteepestDescent, self).__init__(geometry, **kwargs)

        assert(self.alpha > 0), "Alpha should be positive!"

    def optimize(self):
        return self.scale_by_max_step(self.alpha*self.forces[-1])
