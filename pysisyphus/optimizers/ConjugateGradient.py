from typing import Literal

import numpy as np

from pysisyphus.optimizers.BacktrackingOptimizer import BacktrackingOptimizer

# http://ikuz.eu/2015/04/15/the-concept-of-conjugate-gradient-descent-in-python/


class ConjugateGradient(BacktrackingOptimizer):
    def __init__(
        self,
        geometry,
        alpha=0.1,
        formula: Literal["FR", "PR"] = "FR",
        dont_skip=True,
        **kwargs,
    ):
        super().__init__(geometry, alpha=alpha, **kwargs)

        self.formula = formula
        self.dont_skip = dont_skip

    def reset(self):
        super().reset()

        # Check if the number of coordinates changed
        if self.forces[-1].size != self.geometry.coords.size:
            new_forces = self.geometry.forces
            self.forces.append(new_forces)
        self.resetted = True

    def get_beta(self, cur_forces, prev_forces):
        # Fletcher-Reeves formula
        if self.formula == "FR":
            beta = cur_forces.dot(cur_forces) / prev_forces.dot(prev_forces)
        # Polak-Ribiere
        elif self.formula == "PR":
            beta = -cur_forces.dot(prev_forces - cur_forces) / prev_forces.dot(
                prev_forces
            )
            beta_old = cur_forces.dot(cur_forces - prev_forces) / prev_forces.dot(
                prev_forces
            )
            self.log(f"beta_old={beta_old:.4f}, beta={beta:.4f}")
            if beta < 0:
                self.log(f"beta = {beta:.04f} < 0, resetting to 0")
                # beta = 0 basically restarts CG, as no previous step
                # information is mixed into the current step.
                beta = 0
        return beta

    def get_step(
        self, energy, forces, hessian=None, eigvals=None, eigvecs=None, resetted=None
    ):
        if self.cur_cycle > 0:
            self.skip = self.backtrack(forces, self.forces[-2])
        """
        # TODO: do sth. with skip?!
        skip = self.backtrack(new_forces, cur_forces)
        # Imho backtracking gives bad results here, so only use it if
        # explicitly requested (self.dont_skip == False).
        if (not self.dont_skip) and skip:
            self.geometry.coords = cur_coords
            return None
        """

        if not self.resetted and self.cur_cycle > 0:
            beta = self.get_beta(forces, self.forces[-2])
            self.log(f"beta = {beta:.06f}")
            if np.isinf(beta):
                beta = 1.0
            step = forces + beta * self.steps[-1]
        else:
            # Start with steepest descent in the first iteration
            step = forces
            self.resetted = False

        step = self.alpha * step
        step = self.scale_by_max_step(step)
        return step
