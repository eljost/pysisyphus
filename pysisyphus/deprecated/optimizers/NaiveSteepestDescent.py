#!/usr/bin/env python3

import numpy as np

from pysisyphus.helpers import procrustes
from pysisyphus.optimizers.BacktrackingOptimizer import BacktrackingOptimizer

class NaiveSteepestDescent(BacktrackingOptimizer):

    def __init__(self, geometry, **kwargs):
        super(NaiveSteepestDescent, self).__init__(geometry, alpha=0.1, **kwargs)

    def optimize(self):
        if self.is_cos and self.align:
            procrustes(self.geometry)

        self.forces.append(self.geometry.forces)

        step = self.alpha*self.forces[-1]
        step = self.scale_by_max_step(step)
        return step
