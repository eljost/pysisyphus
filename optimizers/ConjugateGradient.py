#!/usr/bin/env python3

import numpy as np

from pysisyphus.optimizers.BacktrackingOptimizer import BacktrackingOptimizer

# http://ikuz.eu/2015/04/15/the-concept-of-conjugate-gradient-descent-in-python/

class ConjugateGradient(BacktrackingOptimizer):

    def __init__(self, geometry, **kwargs):
        self.alpha = -0.05
        super(ConjugateGradient, self).__init__(geometry, **kwargs)

        assert(self.alpha < 0), "Alpha should be negative!"

    def optimize(self):
        if self.cur_cycle > 0:
            beta = (np.vdot(self.forces[-1], self.forces[-1]) /
                    np.vdot(self.forces[-2], self.forces[-2])
            )
            steps = self.forces[-1] + beta*self.steps[-1]
            #print("BETA", beta)
        else:
            steps = self.forces[-1]
        
        steps = self.alpha*steps
        steps = self.scale_by_max_step(steps)
        #new_coords = 

        #    self.skip = self.backtrack()
        return steps
