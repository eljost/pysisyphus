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

    def _get_opt_restart_info(self):
        opt_restart_info = {
            "alpha": self.alpha,
            "cycles_since_backtrack": self.cycles_since_backtrack,
        }
        return opt_restart_info

    def _set_opt_restart_info(self, opt_restart_info):
        self.alpha = opt_restart_info["alpha"]
        self.cycles_since_backtrack = opt_restart_info["cycles_since_backtrack"]

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
