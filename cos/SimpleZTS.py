#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import splprep, splev

from cos.ChainOfStates import ChainOfStates

# [1] http://aip.scitation.org/doi/abs/10.1063/1.2720838

class SimpleZTS(ChainOfStates):

    def __init__(self, images):
        self.alpha = -0.05
        super(SimpleZTS, self).__init__(images)

        self.all_coords = list()
        self.all_forces = list()

    @property
    def forces(self):
        inner_indices = list(range(1, len(self.images)-1))
        [self.images[i].calc_energy_and_forces() for i in inner_indices]
        self.fix_endpoints()
        self._forces  = np.concatenate([image.forces for image in self.images])
        self.all_forces.append(self._forces)
        return self._forces

    def run(self):
        self.all_coords.append(self.coords)
        # [1], Eq. (10)
        steps = self.alpha*self.forces
        self.coords = self.coords + steps

        # Reparametrize mesh
        # To use splprep we have to transpose the coords. 
        transp_coords = np.array([image.coords for image in self.images]).transpose()
        uniform_mesh = np.linspace(0, 1, num=len(self.images))
        tck, u = splprep(transp_coords, s=0)
        new_points = np.array(splev(uniform_mesh, tck))
        new_points = new_points.transpose()
        self.coords = new_points.flatten()
