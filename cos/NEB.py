#!/usr/bin/env python3

import numpy as np

from pysisyphus.cos.ChainOfStates import ChainOfStates

# [1] http://aip.scitation.org/doi/pdf/10.1063/1.1323224
# 
# https://github.com/cstein/neb/blob/master/neb/neb.py


class NEB(ChainOfStates):

    def __init__(self, images):
        super(NEB, self).__init__(images)

        self.perp_forces = list()
        self.par_forces = list()

    @property
    def parallel_forces(self):
        indices = range(len(self.images))
        par_forces = [self.get_parallel_forces(i) for i in indices]
        return np.array(par_forces).flatten()

    def get_parallel_forces(self, i):
        k = 0.1
        if (i is 0) or (i is len(self.images) - 1):
            return k * self.get_tangent(i)

        prev_image = self.images[i-1].coords
        image = self.images[i].coords
        next_image = self.images[i+1].coords
        return (k * (np.linalg.norm(next_image-image) -
                np.linalg.norm(image-prev_image)) *
                self.get_tangent(i)
        )

    @property
    def forces(self):
        indices = range(len(self.images))

        # The calculation of the forces for every image is dependent on
        # the energies of its neighbouring images via the calculation of
        # the tangents.
        if self._forces is None:
            [image.calc_energy_and_forces() for image in self.images]

        perp_forces = np.array(
            [self.get_perpendicular_forces(i) for i in indices]
        )
        par_forces = np.array(
            [self.get_parallel_forces(i) for i in indices]
        )
        total_forces = perp_forces + par_forces
        self._forces = total_forces.flatten()

        return self._forces
