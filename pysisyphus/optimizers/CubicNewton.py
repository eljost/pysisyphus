# [1] https://arxiv.org/abs/2112.02089
#     Regularized Newton Method with Global O(1/k²) Convergence
#     Konstantin Mishchenko

from math import sqrt

import numpy as np

from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer
from pysisyphus.optimizers.exceptions import OptimizationError


class CubicNewton(HessianOptimizer):
    def __init__(self, geometry, **kwargs):
        # Force-disable trust radius update, as this optimizers uses a line search
        kwargs["trust_update"] = False
        super().__init__(geometry, **kwargs)

        self.line_search_cycles = 0

    def optimize(self):
        energy, gradient, hessian, *_ = self.housekeeping()

        if self.cur_cycle == 0:
            # Initial Lipschitz constant estimate; line 2 in algorithm 2 in [1]
            trial_step_length = 0.1
            trial_step = trial_step_length * (-gradient / np.linalg.norm(gradient))
            trial_coords = self.geometry.coords + trial_step
            trial_results = self.geometry.get_energy_and_forces_at(trial_coords)
            trial_gradient = -trial_results["forces"]
            H = (
                np.linalg.norm(trial_gradient - gradient - hessian.dot(trial_step))
                / np.linalg.norm(trial_step) ** 2
            )
        else:
            H = self.H_prev / 4
        self.log(f"Lipschitz constant in cycle {self.cur_cycle}, H={H:.4f}")

        for i in range(self.max_micro_cycles):
            self.line_search_cycles += 1
            H *= 2
            self.log(f"Adaptive Newton line search, cycle {i} using H={H:.4f}")
            lambda_ = sqrt(H * np.linalg.norm(gradient))
            # Instead of solving the linear system we could also use the
            # eigenvectors/-values from housekeeping(). Currently, they are
            # gathered in '*_' and not used.
            trial_step = np.linalg.solve(
                hessian + lambda_ * np.eye(gradient.size), -gradient
            )
            trial_step_norm = np.linalg.norm(trial_step)
            trial_coords = self.geometry.coords + trial_step

            trial_results = self.geometry.get_energy_and_forces_at(trial_coords)
            trial_gradient = -trial_results["forces"]
            trial_energy = trial_results["energy"]

            trial_gradient_small_enough = (
                np.linalg.norm(trial_gradient) <= 2 * lambda_ * trial_step_norm
            )
            sufficient_energy_lowering = (
                trial_energy <= energy - 2 / 3 * lambda_ * trial_step_norm ** 2
            )

            if trial_gradient_small_enough and sufficient_energy_lowering:
                step = trial_step
                break
        else:
            raise OptimizationError("Adaptive Newton line search failed!")

        self.H_prev = H

        return step

    def postprocess_opt(self):
        self.log(f"Line search cycles: {self.line_search_cycles}")
