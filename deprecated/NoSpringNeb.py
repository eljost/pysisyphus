#!/usr/bin/env python3

import numpy as np

from pysisyphus.cos.NEB import NEB

class NoSpringNEB(NEB):

    #def __init__(self, calculator, images):
    def __init__(self, calculator, images):
        super(NoSpringNEB, self).__init__(calculator, images)

    def get_perpendicular_force(self, i):
        return np.array(self.calculator.get_grad(*self.images[i]))

    def get_parallel_force(self, i):
        return (0, 0)
