#!/usr/bin/env python3

import logging

import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer

class BacktrackingOptimizer(Optimizer):

    def __init__(self, geometry, alpha, bt_force=5,
                 dont_skip_after=2, bt_max_scale=4,
                 bt_disable=False, **kwargs):
        # Setting some default values
        self.alpha = alpha
        assert(self.alpha > 0), "Alpha must be positive!"
        self.bt_force = bt_force
        self.dont_skip_after = dont_skip_after
        self.bt_max_scale = bt_max_scale
        self.bt_disable = bt_disable
        assert(self.dont_skip_after >= 1)
        self.cycles_since_backtrack = self.bt_force
        self.scale_factor = 0.5

        super(BacktrackingOptimizer, self).__init__(geometry, **kwargs)

        self.alpha0 = self.alpha
        self.alpha_max = self.bt_max_scale * self.alpha0

        # Keep the skipping history to avoid infinite skipping, e.g. always
        # return skip = False if we already skipped in the last n iterations.
        self.skip_log = list()

    def save_also(self):
        return {
            "cycles_since_backtrack": self.cycles_since_backtrack,
            "alpha": self.alpha,
            "alpha0": self.alpha0,
        }

    def scale_alpha(self, unscaled_steps, alpha):
        # When using an accelerated backtracking optimizer we will vary
        # alpha until a suitable step size is found. If we did a bad step
        # that is discared by self.backtrack(...) we will scale down alpha
        # and try a new step from the last geometry.
        #
        # steps = alpha * direction
        #
        # If alpha is unreasonably big, scaling by a maximum step size will
        # lead to the same scaled step even for different alphas as obtained
        # after several skipping backtracking iterations. So scaling alpha
        # down won't yield an improved step, until alpha becomes very small.
        #
        # If the biggest element in the steps vector for alpha = 1 is e.g.
        # 50 times the maximum allwod step size, halving alpha to 0.5 will
        # still yield a steps vector with a maximum element of 25 times the
        # maximum allowed step size. Both steps vectors would yield the
        # same scaled down steps vector that won't improve our geometry.

        scaled_alphas = alpha * self.scale_factor**np.arange(5)
        scaled_steps = unscaled_steps * scaled_alphas[:,None]
        scaled_steps = np.apply_along_axis(self.scale_by_max_step, 1, scaled_steps)
        max_steps = np.apply_along_axis(np.max, 1, scaled_steps)
        # Check which alpha produces steps below the maximum size.
        below_max_step = np.argmax(max_steps < self.max_step)
        new_alpha = scaled_alphas[below_max_step]
        print(f"First improvement is expected for {new_alpha:.06f}")
        self.alpha = new_alpha
        print(f"got alpha {alpha}, will use new alpha {new_alpha}")

    def reset(self):
        if self.alpha > self.alpha0:
            self.alpha = self.alpha0
            self.log(f"Resetting! Current alpha is {self.alpha}. Lowering "
                     f"it to {self.alpha0}.")

    def backtrack(self, cur_forces, prev_forces, reset_hessian=None):
        """Accelerated backtracking line search."""
        if self.bt_disable:
            return False

        epsilon = 1e-3

        rms = lambda f: np.sqrt(np.mean(np.square(f)))
        cur_rms_force = rms(cur_forces)
        prev_rms_force = rms(prev_forces)

        rms_diff = (
            (cur_rms_force - prev_rms_force) /
            np.abs(cur_rms_force + prev_rms_force)
        )

        # Skip tells us if we overshot
        skip = False

        # When the optimiziation is converging cur_forces will
        # be smaller than prev_forces, so rms_diff will be negative
        # and hence smaller than epsilon, which is a positive number.

        # We went uphill, slow alpha
        self.log(f"Backtracking: rms_diff = {rms_diff:.03f}")
        if rms_diff > epsilon:
            self.log(f"Scaling alpha with {self.scale_factor:.03f}")
            # self.alpha = max(self.alpha0*.5, self.alpha*self.scale_factor)
            self.alpha *= self.scale_factor
            skip = True
            self.cycles_since_backtrack = self.bt_force
        # We continnue going downhill, rms_diff is smaller than epsilon
        else:
            self.cycles_since_backtrack -= 1
            # Check if we didn't accelerate in the previous cycles
            if self.cycles_since_backtrack < 0:
                self.cycles_since_backtrack = self.bt_force
                if self.alpha < self.alpha0:
                    # Reset alpha
                    self.alpha = self.alpha0
                    skip = True
                    self.log(f"Reset alpha to alpha0 = {self.alpha0:.4f}")
                else:
                    # Accelerate alpha
                    self.alpha /= self.scale_factor
                    self.log(f"Scaled alpha to {self.alpha:.4f}")

        # Avoid huge alphas
        if self.alpha > self.alpha_max:
            self.alpha = self.alpha_max
            self.log("Didn't accelerate as alpha would become too large. "
                     f"keeping it at {self.alpha}.")

        # Don't skip if we already skipped the previous iterations to
        # avoid infinite skipping.
        if ((len(self.skip_log) >= self.dont_skip_after)
            and all(self.skip_log[-self.dont_skip_after:])):
            self.log(f"already skipped last {self.dont_skip_after} "
                      "iterations don't skip now.")
            skip = False
            if self.alpha > self.alpha0:
                self.alpha = self.alpha0
                self.log("Resetted alpha to alpha0.")
        self.skip_log.append(skip)
        self.log(f"alpha = {self.alpha:.4f}, skip = {skip}")

        return skip
