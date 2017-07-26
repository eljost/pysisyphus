#!/usr/bin/env python3

from optimizers.Optimizer import Optimizer

class SteepestDescent(Optimizer):

    def __init__(self, geometry):
        super(SteepestDescent, self).__init__(geometry)

    def optimize(self):
        alpha = -0.05
        step = alpha*self.forces[-1]
        new_coords = self.geometry.coords + step
        self.geometry.coords = new_coords
