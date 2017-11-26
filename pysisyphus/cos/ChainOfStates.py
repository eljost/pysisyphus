#!/usr/bin/env python3

import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.xyzloader import make_trj_str
from pysisyphus.constants import BOHR2ANG

# [1] http://dx.doi.org/10.1063/1.1323224

class ChainOfStates:

    def __init__(self, images, parallel=0, fix_ends=False,
                 fix_first=False, fix_last=False):
        assert(len(images) >= 2), "Need at least 2 images!"
        self.images = images
        self.parallel = parallel
        self.fix_first = fix_ends or fix_first
        self.fix_last = fix_ends or fix_last
        self.fix_ends = fix_ends

        self._coords = None
        self._forces = None
        self._energy = None
        # For an image of N atoms coords_lengths will be 3N
        self.coords_length = self.images[0].coords.size
        self.zero_vec = np.zeros(self.coords_length)

    @property
    def moving_indices(self):
        """Returns the indices of the images that aren't fixed and can be
        optimized."""
        indices = range(len(self.images))
        if self.fix_first:
            indices = indices[1:]
        if self.fix_last:
            indices = indices[:-1]
        return indices

    def zero_fixed_vector(self, vector):
        if self.fix_first:
            vector[0] = self.zero_vec
        if self.fix_last:
            vector[-1] = self.zero_vec
        return vector

    def clear(self):
        self._energy = None
        self._forces = None
        self._hessian = None

    def set_vector(self, name, vector, clear=False):
        vec_per_image = vector.reshape(-1, self.coords_length)
        assert(len(self.images) == len(vec_per_image))
        for i in self.moving_indices:
            setattr(self.images[i], name, vec_per_image[i])
        if clear:
            self.clear()

    @property
    def coords(self):
        """Return a flat 1d array containing the coordinates of all images."""
        all_coords = [image.coords for image in self.images]
        self._coords = np.concatenate(all_coords)
        return self._coords

    @coords.setter
    def coords(self, coords):
        """Distribute the flat 1d coords array over all images."""
        self.set_vector("coords", coords, clear=True)

    def set_coords_at(self, i, coords):
        if i in self.moving_indices:
            self.images[i].coords = coords
        # When dealing with a fixed image don't set coords through the
        # property, which would result in resetting the image's caluclated
        # data. Instead assign coords directly.
        else:
            self.images[i]._coords = coords

    @property
    def energy(self):
        self._energy = np.array([image.energy for image in self.images])
        return self._energy

    @energy.setter
    def energy(self, energies):
        """This is needed for some optimizers like CG and BFGS."""
        assert(len(self.images) == len(energies))
        for i in self.moving_indices:
            self.images[i].energy = energies[i]

        self._energy = energies

    @property
    def forces(self):
        forces = [image.forces for image in self.images]
        forces = self.zero_fixed_vector(forces)
        self._forces  = np.concatenate(forces)
        return self._forces

    @forces.setter
    def forces(self, forces):
        self.set_vector("forces", forces)

    @property
    def perpendicular_forces(self):
        indices = range(len(self.images))
        perp_forces = [self.get_perpendicular_forces(i) for i in indices]
        return np.array(perp_forces).flatten()

    @property
    def masses_rep(self):
        return np.array([image.masses_rep for image in self.images]).flatten()

    @property
    def results(self):
        tmp_results = list()
        for image in self.images:
            res = image.results
            res["coords"] = image.coords
            tmp_results.append(res)
        return tmp_results

    def interpolate_between(self, initial_ind, final_ind, image_num):
        # Check for atom ordering
        initial_image = self.images[initial_ind]
        final_image = self.images[final_ind]
        if not initial_image.atoms == final_image.atoms:
            raise Exception("Wrong atom ordering between images. "
                            "Check your input files!")
        initial_coords = initial_image.coords
        final_coords = final_image.coords
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
        """ [1] Equations (8) - (11)"""
        prev_index = max(i - 1, 0)
        next_index = min(i + 1, len(self.images)-1)

        prev_image = self.images[prev_index]
        ith_image = self.images[i]
        next_image = self.images[next_index]

        prev_coords = prev_image.coords
        ith_coords = ith_image.coords
        next_coords = next_image.coords

        tangent_plus = next_coords - ith_coords
        tangent_minus = ith_coords - prev_coords

        # Handle first and last image
        if i is 0:
            return tangent_plus/np.linalg.norm(tangent_plus)
        elif i is (len(self.images) - 1):
            return tangent_minus/np.linalg.norm(tangent_minus)

        prev_energy = prev_image.energy
        ith_energy = ith_image.energy
        next_energy = next_image.energy

        next_energy_diff = abs(next_energy - ith_energy)
        prev_energy_diff = abs(prev_energy - ith_energy)
        delta_energy_max = max(next_energy_diff, prev_energy_diff)
        delta_energy_min = min(next_energy_diff, prev_energy_diff)

        # Uphill
        if next_energy > ith_energy > prev_energy:
            tangent = tangent_plus
        # Downhill
        elif next_energy < ith_energy < prev_energy:
            tangent = tangent_minus
        # Minimum or Maximum
        else:
            if next_energy >= prev_energy:
                tangent = (tangent_plus * delta_energy_max +
                           tangent_minus * delta_energy_min
                )
            # next_energy < prev_energy
            else:
                tangent = (tangent_plus * delta_energy_min +
                           tangent_minus * delta_energy_max
                )

        normalized_tangent = tangent/np.linalg.norm(tangent)
        return normalized_tangent

    def get_tangents(self):
        return np.array([self.get_tangent(i)
                         for i in range(len(self.images))]
        )

    def get_perpendicular_forces(self, i):
        """ [1] Eq. 12"""
        # Our goal in optimizing a ChainOfStates is minimizing the
        # perpendicular force. Alaways return zero perpendicular
        # forces for fixed images, so that they don't interfere
        # with the convergence check.
        if i not in self.moving_indices:
            return self.zero_vec

        forces = self.images[i].forces
        tangent = self.get_tangent(i)
        return forces - (np.dot(forces, tangent)*tangent)

    def as_xyz(self):
        atoms = self.images[0].atoms
        coords_list = [image.coords.reshape((-1,3)) * BOHR2ANG
                       for image in self.images]
        trj_str = make_trj_str(atoms, coords_list)
        return trj_str
