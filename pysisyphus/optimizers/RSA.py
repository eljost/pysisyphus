from math import sqrt

import numpy as np
from scipy.optimize import root_scalar

from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer


class RSA(HessianOptimizer):
    """The Importance of Step Control in Optimization Methods, del Campo, 2009."""

    def optimize(self):
        energy, gradient, H, big_eigvals, big_eigvecs, resetted = self.housekeeping()

        assert big_eigvals.argmin() == 0
        min_eigval = big_eigvals[0]
        pos_definite = min_eigval > 0.0
        gradient_trans = big_eigvecs.T.dot(gradient)
        # This will be also be True when we come close to a minimizer,
        # but then the Hessian will also be positive definite and a
        # simple Newton step will be used.
        hard_case = abs(gradient_trans[0]) <= 1e-6
        self.log(f"Smallest eigenvalue: {min_eigval:.6f}")
        self.log(f"Positive definite Hessian: {pos_definite}")
        self.log(f"Hard case: {hard_case}")

        def get_step(lambda_):
            return -gradient_trans / (big_eigvals + lambda_)

        # Unshifted Newton step
        newton_step = get_step(0.0)
        newton_norm = np.linalg.norm(newton_step)

        # def on_trust_radius(step, thresh=1e-3):
        # return abs(self.trust_radius - np.linalg.norm(step)) <= thresh

        def on_trust_radius_lin(step):
            return 1 / self.trust_radius - 1 / np.linalg.norm(step)

        def finalize_step(lambda_):
            step = get_step(lambda_)
            step = big_eigvecs.dot(step)
            predicted_change = step.dot(gradient) + 0.5 * step.dot(H).dot(step)
            self.predicted_energy_changes.append(predicted_change)
            return step

        # Simplest case. Positive definite Hessian and predicted step is
        # already in trust radius.
        if pos_definite and newton_norm <= self.trust_radius:
            lambda_ = 0.0
            self.log("Using unshifted Newton step.")
            return finalize_step(lambda_)

        # If the Hessian is not positive definite or if the step is too
        # long we have to determine the shift parameter lambda.
        rs_kwargs = {
            "f": lambda lambda_: on_trust_radius_lin(get_step(lambda_)),
            "xtol": 1e-3,
            # Would otherwise be chosen automatically, but we set it
            # here explicitly for verbosity.
            "method": "brentq",
        }

        def root_search(bracket):
            rs_kwargs.update(
                {
                    "bracket": bracket,
                    "x0": bracket[0] + 1e-3,
                }
            )
            res = root_scalar(**rs_kwargs)
            return res

        BRACKET_END = 1e10
        if not hard_case:
            bracket_start = 0.0 if pos_definite else -min_eigval - 1e-3
            bracket = (bracket_start, BRACKET_END)
            res = root_search(bracket)
            assert res.converged
            return finalize_step(res.root)

        # Hard case.
        # First we try the bracket (-b1, ∞)
        bracket = (-min_eigval, BRACKET_END)
        res = root_search(bracket)
        if res.converged:
            return finalize_step(res.root)

        # Now we would try the bracket (-b2, -b1). The resulting step should have
        # a suitable length, but the (shifted) Hessian would have an incorrect
        # eigenvalue spectrum (not positive definite). To solve this we use a
        # different formula to calculate the step.
        without_first = gradient_trans[1:] / (big_eigvals[1:] - min_eigval)
        tau = sqrt(self.trust_radius ** 2 - (without_first ** 2).sum())
        step_trans = [tau] + -without_first.tolist()
        step = big_eigvecs.dot(step_trans)
        predicted_change = step.dot(gradient) + 0.5 * step.dot(H).dot(step)
        self.predicted_energy_changes.append(predicted_change)
        return step
