#/!usr/bin/env python3

# See [1] 10.1063/1.1691018

import numpy as np
from scipy.interpolate import splprep, splev

from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.cos.GrowingChainOfStates import GrowingChainOfStates


class GrowingString(GrowingChainOfStates):

    def __init__(self, images, calc_getter, perp_thresh=0.05,
                 reparam_every=3, **kwargs):
        assert len(images) >= 2, "Need at least 2 images for GrowingString."
        if len(images) > 2:
            images = [images[0], images[-1]]
            print("More than 2 images were supplied! Will only use the "
                  "first and last images to start the GrowingString."
            )
        super().__init__(images, calc_getter, **kwargs)
        self.perp_thresh = perp_thresh
        self.reparam_every = int(reparam_every)
        assert self.reparam_every >= 1

        left_img, right_img = self.images

        self.left_string = [left_img, ]
        self.right_string = [right_img, ]

        self.reparam_in = reparam_every
        self._tangents = None
        self.tangent_list = list()
        self.perp_forces_list = list()
        self.coords_list = list()

        # Adding nodes requires a direction/tangent. As we start
        # from two images we can't do a cubic spline yet, so there are also
        # no splined tangents.
        # For the first two new nodes we use a simple tangent: the
        # unit vector pointing from the left to the right image.
        # As "get_tangent" is reimplemented in this class we call the
        # implementation of the parent ChainOfStates class.
        init_tangent = super().get_tangent(0)
        # Initial distances between left and right image
        Sk, _ = self.arc_dims
        S = Sk / (self.max_nodes+1)
        # Create first two mobile nodes
        left_img, right_img = self.images
        new_left_coords = left_img.coords + S*init_tangent
        new_right_coords = right_img.coords - S*init_tangent
        left_frontier = self.get_new_image(new_left_coords, 1)
        self.left_string.append(left_frontier)
        right_frontier = self.get_new_image(new_right_coords, 2)
        self.right_string.append(right_frontier)

        # Now we have four images and can calculate an initial set of tangents
        # as first derivative of the cubic spline.
        self.set_tangents()
        # The desired spacing of the nodes in the final string on the
        # normalized arclength.
        self.sk = 1 / (self.max_nodes+1)

    @property
    def left_size(self):
        return len(self.left_string)

    @property
    def right_size(self):
        return len(self.right_string)

    @property
    def string_size(self):
        return self.left_size + self.right_size

    @property
    def images_left(self):
        """Returns wether we already created all images."""
        return (self.left_size-1 + self.right_size-1) < self.max_nodes

    @property
    def lf_ind(self):
        """Index of the left frontier node in self.images."""
        return len(self.left_string)-1

    @property
    def rf_ind(self):
        """Index of the right frontier node in self.images."""
        return self.lf_ind+1

    def spline(self):
        reshaped = self.coords.reshape(-1, self.coords_length)
        # To use splprep we have to transpose the coords.
        transp_coords = reshaped.transpose()
        # Spline in batches as scipy can't handle > 11 rows at once
        tcks, us = zip(*[splprep(transp_coords[i:i+9], s=0, k=3)
                         for i in range(0, len(transp_coords), 9)]
        )
        return tcks, us

    def reparam(self, tcks, param_density):
        # Reparametrize mesh
        new_points = np.vstack([splev(param_density, tck) for tck in tcks])
        # Flatten along first dimension.
        new_points = new_points.reshape(-1, len(self.images))
        self.coords = new_points.transpose().flatten()

    def set_tangents(self):
        """Set tangents as given by the first derivative of a cubic spline.

        This method may be considerd hacky as it calculates all
        tangents at once, and not one-by-one as the parent class
        implementation.

        !!! Right now one must not forget to call this method
        after coordinate modification, e.g. after
        reparametrization!
        Otherwise some wrong old tangets are used.
        """

        self.log("Setting tangents using hacky method :)")
        tcks, us = self.spline()
        Sk, cur_mesh = self.arc_dims
        self.log(f"Total arclength Sk={Sk:.4f}")
        tangents = np.vstack([splev(cur_mesh, tck, der=1) for tck in tcks]).T
        norms = np.linalg.norm(tangents, axis=1)
        tangents = tangents / norms[:,None]
        # Tangents of the right string shall point towards the center, so
        # we reverse their orientation.
        tangents[self.rf_ind:] *= -1
        self._tangents = tangents

    def get_tangent(self, i):
        return self._tangents[i]

    @ChainOfStates.forces.getter
    def forces(self):
        if self._forces is None:
            self.calculate_forces()

        indices = range(len(self.images))
        perp_forces = [self.get_perpendicular_forces(i) for i in indices]
        self._forces = np.array(perp_forces).flatten()
        self.perp_forces_list.append(self._forces.copy())
        # TODO: Add climbing modification
        # total_forces = self.set_climbing_forces(total_forces)

        return self._forces

    def reparametrize(self):
        # Non-fully-grown strings are reparametrized every cycle
        can_reparametrize = True
        self.reparam_in -= 1
        # Fully-grown strings are reparametrized only every n-th cycle
        if (self.images_left == 0) and not (self.reparam_in == 0):
            self.log("Skipping reparametrization. Next reparametrization in "
                     f"{self.reparam_in} cycles.")
            self.set_tangents()
            return False

        # Spline displaced coordinates
        tcks, us = self.spline()

        # Calculate the norm of the perpendicular force for every
        # node/image on the string.
        perp_forces  = self.perp_forces_list[-1].reshape(len(self.images), -1)
        perp_norms = np.linalg.norm(perp_forces, axis=1)

        # We can add a new node if the norm of the perpendicular force
        # on the frontier node(s) is below a threshold.
        def converged(i):
            return perp_norms[i] <= self.perp_thresh

        # Check if we can add new nodes
        if self.images_left and converged(self.lf_ind):
            # Insert at the end of the left string, just before the
            # right frontier node.
            new_left_frontier = self.get_new_image(self.dummy_coords, self.rf_ind)
            self.left_string.append(new_left_frontier)
            self.log("Added new left frontier node.")
        if self.images_left and converged(self.rf_ind):
            # Insert at the end of the right string, just before the
            # current right frontier node.
            new_right_frontier = self.get_new_image(self.dummy_coords, self.rf_ind)
            self.right_string.append(new_right_frontier)
            self.log("Added new right frontier node.")

        self.log(f"Current string size: {self.left_size}+{self.right_size}")
        # Reparametrize nodes
        left_inds = np.arange(self.left_size)
        right_inds = np.arange(self.max_nodes+2)[-self.right_size:]
        param_inds = np.concatenate((left_inds, right_inds))
        param_density = self.sk*param_inds
        self.log(f"New param density: " + np.array2string(param_density, precision=2))
        self.reparam(tcks, param_density)

        self.set_tangents()
        self.reparam_in = self.reparam_every

        return True

    def get_additional_print(self):
        size_str = f"{self.left_size}+{self.right_size}"
        if self.images_left == 0:
            size_str = "Full"
        size_info = f"String={size_str}"
        energies = np.array(self.all_energies[-1])
        barrier = (energies.max() - energies[0]) * AU2KJPERMOL
        barrier_info = f"(E_max-E_0)={barrier:.1f} kJ/mol"
        hei_ind = energies.argmax()
        hei_str = f"HEI={hei_ind+1}/{energies.size}"

        tot = f"Grads={self.get_image_calc_counter_sum()}"

        strs = (
            size_info,
            hei_str,
            barrier_info,
        )
        return "\t" + " ".join(strs)
