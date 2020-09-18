import numpy as np

from pysisyphus.irc.IRC import IRC
from pysisyphus.TableFormatter import TableFormatter


class Euler(IRC):

    def __init__(self, geometry, step_length=0.01, **kwargs):
        super(Euler, self).__init__(geometry, step_length, **kwargs)

        step_header = "E/au ΔE(TS)/au |gradient|".split()
        step_fmts = [".6f", ".6f", ".4f"]
        self.step_formatter = TableFormatter(step_header, step_fmts, 10)

    def step(self):
        grad = self.mw_gradient
        grad_norm = np.linalg.norm(grad)

        # Step downhill, against the gradient
        step_direction = -grad / grad_norm
        self.mw_coords += self.step_length*step_direction
