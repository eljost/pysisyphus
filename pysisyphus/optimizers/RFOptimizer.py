#!/usr/bin/env python3

import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer

class RFOptimizer(Optimizer):

    def __init__(self, geometry, **kwargs):
        super().__init__(geometry, **kwargs)

    def prepare_opt(self):
        self.inv_hessian = self.geometry.get_initial_hessian()

    def scale_by_max_step(self, steps):
        return steps

    def optimize(self):
        pass
        #self.forces.append(self.geometry.forces)
        #self.energies.append(self.geometry.energy)

        #return step
