#!/usr/bin/env python3

import numpy as np

from cos.ChainOfStates import ChainOfStates
from optimizer.steepest_descent import steepest_descent

# [1] http://aip.scitation.org/doi/pdf/10.1063/1.1323224
# 
# https://github.com/cstein/neb/blob/master/neb/neb.py


class NEB(ChainOfStates):

    def __init__(self, calculator, images):
        super(NEB, self).__init__(calculator, images)

    def get_tangent(self, i):
        # [1], Eq. (2)
        prev_image = self.images[i-1]
        next_image = self.images[i+1]
        return (next_image-prev_image) / np.linalg.norm(next_image-prev_image)

    def make_tangents(self):
        self.tangents = np.array([self.get_tangent(i) for i
                                  in range(1, len(self.images)-1)])

    def get_perpendicular_force(self, i):
        grad = np.array(self.calculator.get_grad(*self.images[i]))
        tangent = self.get_tangent(i)
        return grad - (np.vdot(grad, tangent)*tangent)

    def get_parallel_force(self, i):
        k = 0.1
        prev_image = self.images[i-1]
        image = self.images[i]
        next_image = self.images[i+1]
        return (k * (np.linalg.norm(next_image-image) -
                np.linalg.norm(image-prev_image)) *
                self.get_tangent(i)
        )

    def take_step(self):
        self.grad_xs, self.grad_ys = self.calculator.get_grad(self.images[:, 0],
                                                              self.images[:, 1])
        
        inner_indices = list(range(1, len(self.images)-1))
        self.make_tangents()
        self.perp_forces = np.array([self.get_perpendicular_force(i) for i in inner_indices])
        self.par_forces = np.array([self.get_parallel_force(i) for i in inner_indices])
        self.total_forces = self.perp_forces + self.par_forces

        self.old_images = self.images.copy()

        new_images = ([steepest_descent(self.images[i], tf)
                               for i, tf in zip(inner_indices, self.total_forces)]
        )
        self.images = np.vstack((self.images[0].copy(), new_images, self.images[-1].copy()))
