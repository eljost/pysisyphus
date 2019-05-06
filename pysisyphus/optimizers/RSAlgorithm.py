#!/usr/bin/env python3

# [1] The Importance of Step Control in Optimization Methods


import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer
from pysisyphus.optimizers.hessian_updates import bfgs_update


class RSAlgorithm(Optimizer):

    def __init__(self, geometry, trust_radius=0.3, **kwargs):
        super().__init__(geometry, **kwargs)

        self.trust_radius = trust_radius

    def prepare_opt(self):
        self.H = self.geometry.get_initial_hessian()

    def optimize(self):
        gradient = self.geometry.gradient
        self.forces.append(-self.geometry.gradient)
        self.energies.append(self.geometry.energy)

        eigvals, eigvecsT = np.linalg.eigh(self.H)

        gradient_transformed = eigvecsT.T.dot(gradient)

        newton_step = -np.linalg.pinv(self.H).dot(gradient)
        newton_norm = np.linalg.norm(newton_step)
        if newton_norm > self.trust_radius:
            return newton_step

        return step
