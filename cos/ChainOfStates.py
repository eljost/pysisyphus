#!/usr/bin/env python3

import numpy as np

from Geometry import Geometry
from qchelper.geometry import make_trj_str

class ChainOfStates:

    def __init__(self, images):
        self.images = images
        self._coords = None
        self._forces = None
        self.coords_length = self.images[0].coords.size

    @property
    def coords(self):
        """Return one big 1d array containing coordinates of all images."""
        all_coords = [image.coords for image in self.images]
        self._coords = np.concatenate(all_coords)
        return self._coords

    @coords.setter
    def coords(self, coords):
        """Distribute the big 1d coords array over all images."""
        coords = coords.reshape(-1, self.coords_length)
        for image, c in zip(self.images, coords):
            image.coords = c

    @property
    def forces(self):
        raise Exception("Not implemented!")

    def interpolate_between(self, initial_ind, final_ind, image_num):
        initial_coords = self.images[initial_ind].coords
        final_coords = self.images[final_ind].coords
        step = (final_coords-initial_coords) / (image_num+1)
        # initial + i*step
        i_array = np.arange(1, image_num)
        atoms = self.images[0].atoms
        new_coords = initial_coords + i_array[:, None]*step
        return [Geometry(atoms, nc) for nc in new_coords]

    def interpolate(self, image_num=5):
        new_images = list()
        # Iterate over image pairs (i, i+1) and interpolate between them
        for i in range(len(self.images)-1):
            interpol_images = self.interpolate_between(i, i+1, image_num)
            new_images.append(self.images[i])
            new_images.extend(interpol_images)
        # As we only added the i-th image and the new images we have to add
        # the last (i+1)-th image at the end.
        new_images.append(self.images[-1])
        self.images = new_images

    def save(self, out_fn):
        atoms = self.images[0].atoms
        coords_list = [image.coords.reshape((-1,3)) for image in self.images]
        trj_str = make_trj_str(atoms, coords_list)
        with open(out_fn, "w") as handle:
            handle.write(trj_str)

    def fix_endpoints(self):
        zero_forces = np.zeros_like(self.images[0].coords)
        self.images[0].forces = zero_forces
        self.images[-1].forces = zero_forces
