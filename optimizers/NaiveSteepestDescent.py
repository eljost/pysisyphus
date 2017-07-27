#!/usr/bin/env python3

from optimizers.Optimizer import Optimizer

class NaiveSteepestDescent(Optimizer):

    def __init__(self, geometry, **kwargs):
        self.alpha = -0.05
        super(NaiveSteepestDescent, self).__init__(geometry, **kwargs)

        assert(self.alpha < 0), "Alpha should be negative!"

    def optimize(self):
        return self.alpha*self.forces[-1]
