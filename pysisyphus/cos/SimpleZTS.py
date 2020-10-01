# [1] http://aip.scitation.org/doi/abs/10.1063/1.2720838

import numpy as np
from scipy.interpolate import splprep, splev

from pysisyphus.cos.ChainOfStates import ChainOfStates

class SimpleZTS(ChainOfStates):

    def __init__(self, images, param="equal", **kwargs):
        self.param = param
        super(SimpleZTS, self).__init__(images, **kwargs)

    def reparametrize(self):
        def weight_function(mean_energies):
            mean_energies = np.abs(mean_energies)
            weights = mean_energies / mean_energies.max()
            weights = np.sqrt(weights)
            return weights

        reshaped = self.coords.reshape(-1, self.coords_length)
        # To use splprep we have to transpose the coords.
        transp_coords = reshaped.transpose()

        # Equal arc length parametrization, default.
        u = None
        # Energy weighted arc length parametrization.
        if self.param == "energy":
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

            u = arc_segments
        # Use chunks of 9 dimension because splprep can handle at max
        # 11 dimensions (as of scipy 1.0.0). Every atoms already yields
        # three dimensions (the coordinates X, Y, Z).
        #
        # For dim <= 11 it would go like this:
        #
        # tck, u = splprep(transp_coords, u=u, s=0)
        # uniform_mesh = np.linspace(0, 1, num=len(self.images))
        # new_points = np.array(splev(uniform_mesh, tck))
        tcks, us = zip(*[splprep(transp_coords[i:i+9], u=u, s=0)
                         for i in range(0, len(transp_coords), 9)]
        )
        uniform_mesh = np.linspace(0, 1, num=len(self.images))
        # Reparametrize mesh
        new_points = np.vstack([splev(uniform_mesh, tck) for tck in tcks])
        # Flatten along first dimension.
        new_points = new_points.reshape(-1, len(self.images))
        self.coords = new_points.transpose().flatten()

        return True
