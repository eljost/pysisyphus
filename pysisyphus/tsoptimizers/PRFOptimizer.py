#!/usr/bin/env python3

# [1] http://aip.scitation.org/doi/10.1063/1.1515483 Optimization review
# [2] https://doi.org/10.1063/1.450914 Trust region method
# [3] 10.1007/978-0-387-40065-5 Numerical optimization
# [4] 10.1007/s00214-016-1847-3 Explorations of some refinements


import numpy as np

from pysisyphus.optimizers.Optimizer import Optimizer


class PRFOptimizer(Optimizer):

    def __init__(self, geometry, root=-1, trust_radius=0.3, **kwargs):
        super().__init__(geometry, **kwargs)

        self.root = root
        self.trust_radius = trust_radius

        self.min_trust_radius = 0.25*self.trust_radius
        self.trust_radius_max = self.trust_radius
        self.predicted_energy_changes = list()
        self.actual_energy_changes = list()

    def optimize(self):
        gradient = self.geometry.gradient
        self.forces.append(-self.geometry.gradient)
        self.energies.append(self.geometry.energy)

        return step
