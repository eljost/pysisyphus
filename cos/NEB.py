#!/usr/bin/env python3

import copy

import numpy as np

from cos.ChainOfStates import ChainOfStates
from optimizer.steepest_descent import steepest_descent

# [1] http://aip.scitation.org/doi/pdf/10.1063/1.1323224
# 
# https://github.com/cstein/neb/blob/master/neb/neb.py


class NEB(ChainOfStates):

    def __init__(self, images):
        super(NEB, self).__init__(images)

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
        
        self.make_tangents()

        self.perp_forces = np.array([self.get_perpendicular_forces(i)
                                     for i in inner_indices])
        self.par_forces = np.array([self.get_parallel_forces(i)
                                    for i in inner_indices])
        self.total_forces = self.perp_forces + self.par_forces

        for ii, tf in zip(inner_indices, self.total_forces):
            self.images[ii].forces = tf
        self._forces  = np.concatenate(self.total_forces)
        return self._forces
