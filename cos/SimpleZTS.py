#!/usr/bin/env python3

import numpy as np
from scipy.interpolate import splprep, splev

from cos.ChainOfStates import ChainOfStates

# [1] http://aip.scitation.org/doi/abs/10.1063/1.2720838

class SimpleZTS(ChainOfStates):

    def __init__(self, images):
        super(SimpleZTS, self).__init__(images)

        inner_indices = list(range(1, len(self.images)-1))

    def reparametrize(self, param="equal"):

        def weight_function(mean_energies):
            mean_energies = np.abs(mean_energies)
            weights = mean_energies / mean_energies.max()
            return weights

        reshaped = self.coords.reshape(-1, self.coords_length)
        # To use splprep we have to transpose the coords.
        transp_coords = reshaped.transpose()

        # Energy weighted arc length parametrization
        if param == "energy":
            energies = [img.energy for img in self.images]
            mean_energies = [(energies[i] + energies[i-1])/2
                             for i in range(1, len(self.images))
            ]
            weights = weight_function(mean_energies)
            arc_segments = [0, ]
            for i in range(1, len(self.images)):
                coord_diff = np.linalg.norm(
                    self.images[i].coords - self.images[i-1].coords
                )
                next_segment = arc_segments[-1] + weights[i-1]*coord_diff
                arc_segments.append(next_segment)
            arc_segments = np.array(arc_segments)
            arc_segments /= arc_segments.max()

            tck, u = splprep(transp_coords, u=arc_segments, s=0)
        # Equal arc length parametrization
        elif param == "equal":
            tck, u = splprep(transp_coords, s=0)

        # Reparametrize mesh
        uniform_mesh = np.linspace(0, 1, num=len(self.images))
        new_points = np.array(splev(uniform_mesh, tck))
        self.coords = new_points.transpose().flatten()
