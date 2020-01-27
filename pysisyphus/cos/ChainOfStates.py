#!/usr/bin/env python3

from copy import copy
import logging
from multiprocessing import Pool

from distributed import Client
import numpy as np
from scipy.interpolate import interp1d, splprep, splev

from pysisyphus.constants import BOHR2ANG
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import get_coords_diffs
from pysisyphus.xyzloader import make_trj_str


# [1] http://dx.doi.org/10.1063/1.1323224


class ChainOfStates:
    logger = logging.getLogger("cos")
    valid_coord_types = "cart dlc".split()

    def __init__(self, images, parallel=0, fix_ends=True,
                 fix_first=False, fix_last=False,
                 climb=False, climb_rms=5e-3,
                 scheduler=None):

        assert(len(images) >= 2), "Need at least 2 images!"
        self.images = list(images)
        self.parallel = parallel
        self.fix_first = fix_ends or fix_first
        self.fix_last = fix_ends or fix_last
        self.fix_ends = fix_ends
        self.climb = climb
        self.climb_rms = climb_rms
        self.scheduler = scheduler

        self._coords = None
        self._forces = None
        self._energy = None
        self.counter = 0
        self.coords_length = self.images[0].coords.size
        self.cart_coords_length = self.images[0].cart_coords.size
        self.zero_vec = np.zeros(self.coords_length)

        self.coords_list = list()
        self.forces_list = list()
        self.all_energies = list()
        self.all_true_forces = list()

        # Start climbing immediateley with climb_rms == -1
        self.started_climbing = self.climb_rms == -1
        if self.started_climbing:
            self.log("Will start climbing immediately.")

        img0 = self.images[0]
        self.image_atoms = copy(img0.atoms)
        self.coord_type = img0.coord_type
        assert self.coord_type in self.valid_coord_types, \
                f"Invalid coord_type! Supported types are: {self.valid_coord_types}"
        assert all([img.coord_type == self.coord_type for img in self.images]), \
                "coord_type of images differ!"
        try:
            self.prim_indices = img0.internal.prim_indices
        except AttributeError:
            self.prim_indices = None

    def log(self, message):
        self.logger.debug(f"Counter {self.counter+1:03d}, {message}")

    def get_fixed_indices(self):
        fixed = list()
        if self.fix_first:
            fixed.append(0)
        if self.fix_last:
            fixed.append(len(self.images)-1)
        return fixed

    @property
    def moving_indices(self):
        """Returns the indices of the images that aren't fixed and can be
        optimized."""
        fixed = self.get_fixed_indices()
        return [i for i in range(len(self.images))
                if i not in fixed]

    @property
    def last_index(self):
        return len(self.images) - 1

    @property
    def moving_images(self):
        return [self.images[i] for i in self.moving_indices]

    def zero_fixed_vector(self, vector):
        fixed = self.get_fixed_indices()
        for i in fixed:
            vector[i] = self.zero_vec
        return vector

    def clear(self):
        self._energy = None
        self._forces = None
        self._hessian = None
        try:
            self._tangents = None
        except AttributeError:
            # TODO: move this to another logging level?!
            self.log("There are no tangents to reset.")

    @property
    def atoms(self):
        atoms_ = self.images[0].atoms
        return len(self.images) * atoms_

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

    @property
    def cart_coords(self):
        """Return a flat 1d array containing the cartesian coordinates of all
        images."""
        return np.concatenate([image.cart_coords for image in self.images])

    @property
    def coords3d(self):
        assert self.images[0].coord_type == "cart"
        return self.coords.reshape(-1, 3)

    def set_coords_at(self, i, coords):
        """Called from helpers.procrustes with cartesian coordinates.
        Then tries to set cartesian coordinate as self.images[i].coords
        which will raise an error when coord_type != "cart".
        """
        assert self.images[i].coord_type == "cart", \
            "ChainOfStates.set_coords_at() has to be reworked to support " \
            "internal coordiantes. Try to set 'align: False' in the 'opt' " \
            "section of the .yaml input file."
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

    def par_image_calc(self, image):
        image.calc_energy_and_forces()
        return image

    def set_images(self, indices, images):
        for ind, image in zip(indices, images):
            self.images[ind] = image

    def calculate_forces(self):
        # Determine the number of images for which we have to do calculations.
        # There may also be calculations for fixed images, as they need an
        # energy value. But every fixed image only needs a calculation once.
        images_to_calculate = self.moving_images
        image_indices = self.moving_indices
        if self.fix_first and (self.images[0]._energy is None):
            images_to_calculate = [self.images[0]] + images_to_calculate
            image_indices = [0] + list(image_indices)
        if self.fix_last and (self.images[-1]._energy is None):
            images_to_calculate = images_to_calculate + [self.images[-1]]
            image_indices = list(image_indices) + [-1]
        assert len(images_to_calculate) <= len(self.images)

        # Parallel calculation with dask
        if self.scheduler:
            client = self.get_dask_client()
            self.log(client)
            image_futures = client.map(self.par_image_calc, images_to_calculate)
            self.set_images(image_indices, client.gather(image_futures))
        # Parallel calculation with multiprocessing
        elif self.parallel > 0:
            with Pool(processes=self.parallel) as pool:
                par_images = pool.map(self.par_image_calc, images_to_calculate)
                self.set_images(image_indices, par_images)
        # Serial calculation
        else:
            [image.calc_energy_and_forces() for image in images_to_calculate]
        self.set_zero_forces_for_fixed_images()
        self.counter += 1

        energies = [image.energy for image in self.images]
        forces = [image.forces for image in self.images]
        self.all_energies.append(energies)
        self.all_true_forces.append(forces)

    @property
    def forces(self):
        self.set_zero_forces_for_fixed_images()
        forces = [image.forces for image in self.images]
        self._forces  = np.concatenate(forces)
        self.counter += 1
        return self._forces

    @forces.setter
    def forces(self, forces):
        self.set_vector("forces", forces)

    @property
    def perpendicular_forces(self):
        indices = range(len(self.images))
        perp_forces = [self.get_perpendicular_forces(i) for i in indices]
        return np.array(perp_forces).flatten()

    def get_perpendicular_forces(self, i):
        """ [1] Eq. 12"""
        # Our goal in optimizing a ChainOfStates is minimizing the
        # perpendicular force. Always return zero perpendicular
        # forces for fixed images, so that they don't interfere
        # with the convergence check.
        if i not in self.moving_indices:
            return self.zero_vec

        forces = self.images[i].forces
        tangent = self.get_tangent(i)
        perp_forces = forces - forces.dot(tangent)*tangent
        return perp_forces

    @property
    def gradient(self):
        return -self.forces

    @gradient.setter
    def gradient(self, gradient):
        self.forces = -gradient

    @property
    def masses_rep(self):
        return np.array([image.masses_rep for image in self.images]).flatten()

    @property
    def results(self):
        tmp_results = list()
        for image in self.images:
            res = image.results
            res["coords"] = image.coords
            res["cart_coords"] = image.cart_coords
            tmp_results.append(res)
        return tmp_results

    def set_zero_forces_for_fixed_images(self):
        """This is always done in cartesian coordinates, independent
        of the actual coord_type of the images as setting forces only
        work with cartesian forces."""
        zero_forces = np.zeros_like(self.images[0].cart_coords)
        if self.fix_first:
            self.images[0].cart_forces = zero_forces
            self.log("Zeroed forces on fixed first image.")
        if self.fix_last:
            self.images[-1].cart_forces = zero_forces
            self.log("Zeroed forces on fixed last image.")

    def get_tangent(self, i, kind="upwinding"):
        """ [1] Equations (8) - (11)"""

        tangent_kinds = "upwinding simple bisect".split()
        assert kind in tangent_kinds, "Invalid tangent kind! Valid " \
            f"tangent kinds are {self.tangent_kinds}"
        prev_index = max(i - 1, 0)
        next_index = min(i + 1, len(self.images)-1)

        prev_image = self.images[prev_index]
        ith_image = self.images[i]
        next_image = self.images[next_index]

        # If (i == 0) or (i == len(self.images)-1) then one
        # of this tangents is zero.
        tangent_plus = next_image - ith_image
        tangent_minus = ith_image - prev_image

        # Handle first and last image
        if i == 0:
            return tangent_plus/np.linalg.norm(tangent_plus)
        elif i == (len(self.images) - 1):
            return tangent_minus/np.linalg.norm(tangent_minus)

        # [1], Eq. (1)
        if kind == "simple":
            tangent = next_image - prev_image
            tangent /= np.linalg.norm(tangent)
            return tangent
        # [1], Eq. (2)
        elif kind == "bisect":
            first_term = tangent_minus / np.linalg.norm(tangent_minus)
            sec_term = tangent_plus / np.linalg.norm(tangent_plus)
            tangent = first_term + sec_term
            tangent /= np.linalg.norm(tangent)
            return tangent

        # Upwinding tangent from now on
        # [1], Eq. (8) and so on
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

    def as_xyz(self, comments=None):
        return "\n".join([image.as_xyz() for image in self.images])

    def get_dask_client(self):
        return Client(self.scheduler, pure=False, silence_logs=False)

    def get_hei_index(self, energies=None):
        """Return index of highest energy image."""
        if energies is None:
            energies = [image.energy for image in self.images]
        return np.argmax(energies)

    def prepare_opt_cycle(self, last_coords, last_energies, last_forces):
        """Implements additional logic in preparation of the next
        optimization cycle.

        Should be called by the optimizer at the beginning of a new
        optimization cycle. Can be used to implement additional logic
        as needed for AdaptiveNEB etc.
        """
        self.coords_list.append(last_coords)
        self.forces_list.append(last_forces)

        # Return False if we don't want to climb or are already
        # climbing.
        already_climbing = self.started_climbing
        if self.climb and not already_climbing:
            self.started_climbing = self.check_for_climbing_start()
            if self.started_climbing:
                self.log("Starting to climb in next iteration.")
                print("Starting to climb in next iteration.")
        return not already_climbing and self.started_climbing

    def rms(self, arr):
        """Root mean square

        Returns the root mean square of the given array.

        Parameters
        ----------
        arr : iterable of numbers

        Returns
        -------
        rms : float
            Root mean square of the given array.
        """
        return np.sqrt(np.mean(np.square(arr)))

    def check_for_climbing_start(self):
        # Only initiate climbing on a sufficiently converged MEP.
        # This can be determined from a supplied threshold for the
        # RMS force (rms_force) or from a multiple of the
        # RMS force convergence threshold (rms_multiple, default).
        rms_forces = self.rms(self.forces_list[-1])
        return rms_forces <= self.climb_rms

    def get_climbing_indices(self):
        # Index of the highest energy image (HEI)
        hei_index = self.get_hei_index()

        move_inds = self.moving_indices
        # Don't climb it not yet enabled or requested.
        if not (self.climb and self.started_climbing):
            climb_indices = tuple()
        # We can do two climbing (C2) neb if the highest energy image (HEI)
        # is in moving_indices but not the first or last item in this list.
        elif hei_index in move_inds[1:-1]:
            climb_indices = (hei_index-1, hei_index+1)
        # Do one image climbing (C1) neb if the HEI is the first or last
        # item in moving_indices.
        elif (hei_index == 1) or (hei_index == move_inds[-1]):
            climb_indices = (hei_index, )
        # Don't climb when the HEI is the first or last image of the whole
        # NEB.
        else:
            climb_indices = tuple()
            self.log("Want to climb but can't. HEI is first or last image!")
        self.log(f"Climbing indices are {climb_indices}.")
        return climb_indices

    def get_climbing_forces(self, ind):
        climbing_image = self.images[ind]
        ci_forces = climbing_image.forces
        tangent = self.get_tangent(ind)
        climbing_forces = ci_forces - 2*ci_forces.dot(tangent)*tangent

        return climbing_forces, climbing_image.energy

    def set_climbing_forces(self, forces):
        # Avoids calling the other methods with their logging output etc.
        if not self.started_climbing:
            return forces

        for i in self.get_climbing_indices():
            climb_forces, climb_en = self.get_climbing_forces(i)
            forces[i] = climb_forces
            self.log(f"Climbing with image {i}, E = {climb_en:.6f} au")
        return forces

    def get_splined_hei(self):
        # Interpolate energies
        cart_coords = np.array([image.cart_coords for image in self.images])
        coord_diffs = get_coords_diffs(cart_coords)
        energies_spline = interp1d(coord_diffs, self.energy, kind="cubic")
        x_fine = np.linspace(0, 1, 500)
        energies_fine = energies_spline(x_fine)
        # Determine index that yields the highest energy
        hei_ind = energies_fine.argmax()
        hei_x = x_fine[hei_ind]
        hei_energy = energies_fine[hei_ind]

        reshaped = cart_coords.reshape(-1, self.cart_coords_length)
        # To use splprep we have to transpose the coords.
        transp_coords = reshaped.transpose()
        tcks, us = zip(*[splprep(transp_coords[i:i+9], s=0, k=3)
                         for i in range(0, len(transp_coords), 9)]
        )

        # Reparametrize mesh
        hei_coords = np.vstack([splev([hei_x, ], tck) for tck in tcks])
        hei_coords = hei_coords.flatten()

        hei_tangent = np.vstack([splev([hei_x, ], tck, der=1) for tck in tcks]).T
        hei_tangent = hei_tangent.flatten()
        hei_tangent /= np.linalg.norm(hei_tangent)
        return hei_coords, hei_energy, hei_tangent

    def get_image_calc_counter_sum(self):
        return sum([image.calculator.calc_counter for image in self.images])

    def __str__(self):
        return self.__class__.__name__
