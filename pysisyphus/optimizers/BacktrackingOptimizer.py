#!/usr/bin/env python3

import logging

import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer

class BacktrackingOptimizer(Optimizer):

    def __init__(self, geometry, alpha, dont_skip_after=2, **kwargs):
        # Setting some default values
        self.alpha = alpha
        assert(self.alpha > 0), "Alpha should be positive!"
        self.force_backtrack_in = 5
        self.dont_skip_after = dont_skip_after
        assert(self.dont_skip_after >= 1)
        self.cycles_since_backtrack = self.force_backtrack_in

        super(BacktrackingOptimizer, self).__init__(geometry, **kwargs)

        self.alpha0 = self.alpha

        # Keep the skipping history to avoid infinite skipping, e.g. always
        # return skip = False if we already skipped in the last n iterations.
        self.skip_log = list()

    def backtrack(self, cur_forces, prev_forces, reset_hessian=None):
        """Accelerated backtracking line search."""

        if self.started_climbing:
            self.log("backtracking disabled when climbing")
            self.alpha = self.alpha0
            if reset_hessian:
                self.reset_hessian()
            return False

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
        self.log(f"backtracking: rms_diff = {rms_diff:.04f}")
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
                    self.log(f"reset alpha to alpha0 = {self.alpha0}")
                else:
                    # Accelerate alpha
                    self.alpha /= scale_factor
                    self.log("scaled alpha")
        self.log(f"alpha = {self.alpha}, skip = {skip}")

        # Don't skip if we already skipped the previous iterations.
        if all(self.skip_log[-self.dont_skip_after:]):
            self.log(f"already skipped last {self.dont_skip_after} "
                      "iterations don't skip now.")
            skip = False
            if self.alpha > self.alpha0:
                self.alpha = self.alpha0
                self.log("resetted alpha to alpha0.")
        self.skip_log.append(skip)

        return skip
