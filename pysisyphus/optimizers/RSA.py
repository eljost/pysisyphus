import numpy as np
from scipy.optimize import root_scalar

from pysisyphus.optimizers.HessianOptimizer import HessianOptimizer


class RSA(HessianOptimizer):

    def optimize(self):
        energy, gradient, H, big_eigvals, big_eigvecs = self.housekeeping()

        min_eigval = big_eigvals[0]
        pos_definite = min_eigval > 0.
        gradient_trans = big_eigvecs.T.dot(gradient)

        def get_step(lambda_):
            return -gradient_trans /(big_eigvals + lambda_)

        # Unshifted Newton step
        newton_step = get_step(0.)
        newton_norm = np.linalg.norm(newton_step)

        def on_trust_radius(step, thresh=1e-3):
            return abs(self.trust_radius - np.linalg.norm(step)) <= thresh

        def on_trust_radius_lin(step, thresh=1e-3):
            return 1/self.trust_radius - 1/np.linalg.norm(step)

        def finalize_step(lambda_):
            step = get_step(lambda_)
            step = big_eigvecs.dot(step)
            predicted_change = step.dot(gradient) + 0.5 * step.dot(H).dot(step)
            self.predicted_energy_changes.append(predicted_change)
            return step

        # Simplest case. Positive definite Hessian and predicted step is
        # already in trust radius.
        if pos_definite and newton_norm <= self.trust_radius:
            lambda_ = 0.
            return finalize_step(lambda_)

        try:
            bracket_start = 0. if pos_definite else -min_eigval
            rs_kwargs = {
                "f": lambda lambda_: on_trust_radius_lin(get_step(lambda_)),
                "bracket": (bracket_start, 1e10),
                "x0": bracket_start + 1e-3,
                "xtol": 1e-3,
            }
            res = root_scalar(**rs_kwargs)
            assert res.flag == "converged"
            lambda_ = res.root
        except AssertionError as err:
            raise err

        return finalize_step(lambda_)

