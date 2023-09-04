import numpy as np

from pysisyphus.irc.IRC import IRC


class Euler(IRC):
    def __init__(self, geometry, step_length=0.01, **kwargs):
        super().__init__(geometry, step_length, **kwargs)

    def step(self):
        grad = self.mw_gradient
        grad_norm = np.linalg.norm(grad)

        # Step downhill, against the gradient
        step_direction = -grad / grad_norm
        self.mw_coords += self.step_length * step_direction
