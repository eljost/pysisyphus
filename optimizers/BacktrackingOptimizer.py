#!/usr/bin/env python3

import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer

class BacktrackingOptimizer(Optimizer):

    def __init__(self, geometry, **kwargs):
        # Setting some default values
        self.force_backtrack_in = 3
        self.cycles_since_backtrack = self.force_backtrack_in

        super(BacktrackingOptimizer, self).__init__(geometry, **kwargs)


    def backtrack(self):
        """Accelerated backtracking line search."""
        epsilon = 1e-3
        alpha0 = -0.05
        scale_factor = 0.5

        prev_rms_force, cur_rms_force = self.rms_forces[-2:]
        # chk
        rms_diff = (
            (cur_rms_force - prev_rms_force) /
            np.abs(cur_rms_force+prev_rms_force)
        )
        skip = False

        # Slow alpha
        if rms_diff > epsilon:
            self.alpha *= scale_factor
            skip = True
            self.cycles_since_backtrack = self.force_backtrack_in
        else:
            self.cycles_since_backtrack -= 1
            #print("cycles_since_backtrack", self.cycles_since_backtrack)
            if self.cycles_since_backtrack < 0:
                self.cycles_since_backtrack = self.force_backtrack_in
                if self.alpha > alpha0:
                    # Reset alpha
                    alpha = alpha0
                    skip = True
                else:
                    # Accelerate alpha
                    self.alpha /= scale_factor
        return skip
