#!/usr/bin/env python3

import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.xyzloader import make_trj_str

class ChainOfStates:

    def __init__(self, images):
        self.images = images

        self._coords = None
        self._forces = None
        self.coords_length = self.images[0].coords.size

    @property
    def coords(self):
        """Return a flat 1d array containing the coordinates of all images."""
        all_coords = [image.coords for image in self.images]
        self._coords = np.concatenate(all_coords)
        return self._coords

    @coords.setter
    def coords(self, coords):
        """Distribute the flat 1d coords array over all images."""
        coords = coords.reshape(-1, self.coords_length)
        for image, c in zip(self.images, coords):
            image.coords = c

    @property
    def energy(self):
        self._energy = [image.energy for image in self.images]
        return self._energy

    @property
    def forces(self):
        forces = [image.forces for image in self.images]
        self._forces  = np.concatenate(forces)
        return self._forces

    @forces.setter
    def forces(self, forces):
        forces = forces.reshape(-1, self.coords_length)
        for image, f in zip(self.images, forces):
            image.forces = f

    @property
    def perpendicular_forces(self):
        indices = range(len(self.images))
        perp_forces = [self.get_perpendicular_forces(i) for i in indices]
        return np.array(perp_forces).flatten()

    def interpolate_between(self, initial_ind, final_ind, image_num):
        initial_coords = self.images[initial_ind].coords
        final_coords = self.images[final_ind].coords
        step = (final_coords-initial_coords) / (image_num+1)
        # initial + i*step
        i_array = np.arange(1, image_num+1)
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

    def fix_ends(self):
        zero_forces = np.zeros_like(self.images[0].coords)
        self.images[0].forces = zero_forces
        self.images[-1].forces = zero_forces

    def get_tangent(self, i):
        # Use a one-sided difference for the first and last image
        if i == 0:
            prev_index = i
            next_index = 1
        elif i == (len(self.images) - 1):
            prev_index = i - 1
            next_index = len(self.images) - 1
        # If i is an inner index use the image before and after i
        else:
            prev_index = i - 1
            next_index = i + 1
        # [1], Eq. (2)
        prev_image = self.images[prev_index].coords
        next_image = self.images[next_index].coords
        return (next_image-prev_image) / np.linalg.norm(next_image-prev_image)

    def get_perpendicular_forces(self, i):
        forces = self.images[i].forces
        tangent = self.get_tangent(i)
        return forces - (np.vdot(forces, tangent)*tangent)

    def save(self, out_fn):
        atoms = self.images[0].atoms
        coords_list = [image.coords.reshape((-1,3)) for image in self.images]
        trj_str = make_trj_str(atoms, coords_list)
        with open(out_fn, "w") as handle:
            handle.write(trj_str)
