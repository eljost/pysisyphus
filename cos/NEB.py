#!/usr/bin/env python3

import numpy as np

from pysisyphus.cos.ChainOfStates import ChainOfStates

# [1] http://aip.scitation.org/doi/pdf/10.1063/1.1323224
# 
# https://github.com/cstein/neb/blob/master/neb/neb.py


class NEB(ChainOfStates):

    def __init__(self, images, k=0.01):
        super(NEB, self).__init__(images)

        self.k = k

        self.perp_forces = list()
        self.par_forces = list()

    @property
    def parallel_forces(self):
        indices = range(len(self.images))
        par_forces = [self.get_parallel_forces(i) for i in indices]
        return np.array(par_forces).flatten()

    def get_parallel_forces(self, i):
        if (i is 0) or (i is len(self.images) - 1):
            return self.k * self.get_tangent(i)

        prev_coords = self.images[i-1].coords
        ith_coords = self.images[i].coords
        next_coords = self.images[i+1].coords
        return (self.k
                * (np.linalg.norm(next_coords-ith_coords)
                -  np.linalg.norm(ith_coords-prev_coords)
                ) * self.get_tangent(i)
        )

    @property
    def forces(self):
        indices = range(len(self.images))

        # The calculation of the forces for every image is dependent on
        # the energies of its neighbouring images via the calculation of
        # the tangents.
        if self._forces is None:
            [image.calc_energy_and_forces() for image in self.images]

        par_forces = np.array(
            [self.get_parallel_forces(i) for i in indices]
        )
        perp_forces = np.array(
            [self.get_perpendicular_forces(i) for i in indices]
        )
        total_forces = par_forces + perp_forces
        self._forces = total_forces.flatten()

        return self._forces
