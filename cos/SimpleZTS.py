#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import splprep, splev

from cos.ChainOfStates import ChainOfStates

# [1] http://aip.scitation.org/doi/abs/10.1063/1.2720838

class SimpleZTS(ChainOfStates):

    def __init__(self, images):
        super(SimpleZTS, self).__init__(images)

    @property
    def forces(self):
        inner_indices = list(range(1, len(self.images)-1))
        [self.images[i].calc_energy_and_forces() for i in inner_indices]
        self.fix_endpoints()
        self._forces  = np.concatenate([image.forces for image in self.images])
        return self._forces

    def reparametrize(self, coords):
        reshaped = coords.reshape(-1, self.coords_length)
        # To use splprep we have to transpose the coords.
        transp_coords = reshaped.transpose()
        uniform_mesh = np.linspace(0, 1, num=len(self.images))
        tck, u = splprep(transp_coords, s=0)
        # Reparametrize mesh
        new_points = np.array(splev(uniform_mesh, tck))
        return new_points.transpose().flatten()
