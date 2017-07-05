#!/usr/bin/env python3

import numpy as np

class ChainOfStates:

    def __init__(self, calculator, images):
        self.calculator = calculator
        self.images = np.array(images)

    def interpolate_images(self, image_num=10):
        initial = self.images[0]
        final = self.images[-1]
        step = (final-initial) / (image_num+1)
        # initial + i*step
        i_array = np.arange(image_num+2)
        self.images = initial + i_array[:, None]*step
