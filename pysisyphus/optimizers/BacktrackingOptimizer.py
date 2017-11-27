#!/usr/bin/env python3

import logging

import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer

class BacktrackingOptimizer(Optimizer):

    def __init__(self, geometry, alpha, **kwargs):
        # Setting some default values
        self.alpha = alpha
        self.force_backtrack_in = 5
        self.cycles_since_backtrack = self.force_backtrack_in

        super(BacktrackingOptimizer, self).__init__(geometry, **kwargs)

        self.alpha0 = self.alpha
        assert(self.alpha > 0), "Alpha should be positive!"


    def backtrack(self, cur_forces, prev_forces):
        """Accelerated backtracking line search."""
        epsilon = 1e-3
        scale_factor = 0.5

        rms = lambda f: np.sqrt(np.mean(np.square(f)))
        cur_rms_force = rms(cur_forces)
        prev_rms_force = rms(prev_forces)

        rms_diff = (
            (cur_rms_force - prev_rms_force) /
            np.abs(cur_rms_force+prev_rms_force)
        )

        # Skip tells us if we overshot
        skip = False

        # When the optimiziation is converging cur_forces will
        # be smaller than prev_forces, so rms_diff will be negative
        # and hence always smaller than epsilon.

        # Slow alpha if we go uphill.
        if rms_diff > epsilon:
            self.alpha *= scale_factor
            skip = True
            self.cycles_since_backtrack = self.force_backtrack_in
        else:
            self.cycles_since_backtrack -= 1
            if self.cycles_since_backtrack < 0:
                self.cycles_since_backtrack = self.force_backtrack_in
                if self.alpha < self.alpha0:
                    # Reset alpha
                    self.alpha = self.alpha0
                    skip = True
                else:
                    # Accelerate alpha
                    self.alpha /= scale_factor
        self.log(f"alpha = {self.alpha}, skip = {skip}")
        return skip
