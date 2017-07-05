#!/usr/bin/env python3

import numpy as np

from cos.ChainOfStates import ChainOfStates

class SimpleZTS(ChainOfStates):

    def __init__(self, calculator, images):
        super(SimpleZTS, self).__init__(calculator, images)

    def take_step(self):
        self.grad_xs, self.grad_ys = self.calculator.get_grad(self.images[:, 0],
                                                              self.images[:, 1])
        self.grad = np.stack((self.grad_xs, self.grad_ys), axis=-1)
        print("grad", self.grad.shape, self.grad)
        print(self.grad_xs, self.grad_xs.shape)
        self.old_images = self.images.copy()
        # Evolution
        # Forward euler, 0.05 step
        step = 0.05
        self.images = self.images - step * self.grad
