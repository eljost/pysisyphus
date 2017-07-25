#!/usr/bin/env python3

import numpy as np

from Geometry import Geometry
from qchelper.geometry import make_trj_str

class ChainOfStates:

    def __init__(self, images):
        self.images = images
        self._coords = None
        self._forces = None
        self.coord_length = self.images[0].coords.size

    @property
    def coords(self):
        inner_coords = [image.coords for image in self.images[1:-1]]
        self._coords = np.concatenate(inner_coords)
        return self._coords

    @coords.setter
    def coords(self, coords):
        coords = coords.reshape(-1, self.coord_length)
        for inner_image, c in zip(self.images[1:-1], coords):
            inner_image.coords = c

    def interpolate_images(self, image_num=10):
        initial = self.images[0].coords
        final = self.images[-1].coords
        step = (final-initial) / (image_num+1)
        # initial + i*step
        i_array = np.arange(image_num+2)
        atoms = self.images[0].atoms
        new_coords = initial + i_array[:, None]*step
        self.images = [Geometry(atoms, nc) for nc in new_coords]

    def save(self, out_fn):
        atoms = self.images[0].atoms
        coords_list = [image.coords.reshape((-1,3)) for image in self.images]
        trj_str = make_trj_str(atoms, coords_list)
        with open(out_fn, "w") as handle:
            handle.write(trj_str)
