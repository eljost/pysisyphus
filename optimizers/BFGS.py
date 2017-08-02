#!/usr/bin/env python3

import numpy as np

from pysisyphus.optimizers.BacktrackingOptimizer import BacktrackingOptimizer

class BFGS(BacktrackingOptimizer):

    def __init__(self, geometry, **kwargs):
        super(BFGS, self).__init__(geometry, **kwargs)

        self.inv_hessian = np.eye(self.geometry.coords.size)

    def optimize(self):
        if len(self.forces) is 1:
            # Cal
            print("first iteration!")
        import sys; sys.exit()
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
