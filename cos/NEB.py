#!/usr/bin/env python3

import numpy as np

from cos.ChainOfStates import ChainOfStates

# [1] http://aip.scitation.org/doi/pdf/10.1063/1.1323224
# 
# https://github.com/cstein/neb/blob/master/neb/neb.py


class NEB(ChainOfStates):

    def __init__(self, images):
        super(NEB, self).__init__(images)

        self.perp_forces = list()
        self.par_forces = list()

    def get_tangent(self, i):
        # [1], Eq. (2)
        prev_image = self.images[i-1].coords
        next_image = self.images[i+1].coords
        return (next_image-prev_image) / np.linalg.norm(next_image-prev_image)

    def make_tangents(self):
        self.tangents = np.array([self.get_tangent(i) for i
                                  in range(1, len(self.images)-1)])

    def get_perpendicular_forces(self, i):
        forces = self.images[i].forces
        tangent = self.get_tangent(i)
        return forces - (np.vdot(forces, tangent)*tangent)

    def get_parallel_forces(self, i):
        k = 0.1
        prev_image = self.images[i-1].coords
        image = self.images[i].coords
        next_image = self.images[i+1].coords
        return (k * (np.linalg.norm(next_image-image) -
                np.linalg.norm(image-prev_image)) *
                self.get_tangent(i)
        )

    @property
    def forces(self):
        inner_indices = list(range(1, len(self.images)-1))
        [self.images[i].calc_energy_and_forces() for i in inner_indices]
        # Fix the educt and product
        zero_forces = np.zeros_like(self.images[0].coords)
        self.images[0].forces = zero_forces
        self.images[-1].forces = zero_forces
        
        self.make_tangents()

        perp_forces = np.array([self.get_perpendicular_forces(i)
                                for i in inner_indices])
        par_forces = np.array([self.get_parallel_forces(i)
                               for i in inner_indices])
        total_forces = perp_forces + par_forces

        # Store the calculated forces so they can be used later on, e.g.
        # for plotting. total_forces is saved in self._forces.
        self.perp_forces.append(perp_forces)
        self.par_forces.append(par_forces)

        for ii, tf in zip(inner_indices, total_forces):
            self.images[ii].forces = tf
        # Here we also use the zero forces of educt and product
        self._forces  = np.concatenate([image.forces for image in self.images])
        return self._forces
