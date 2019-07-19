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

        left_img, right_img = self.images
        # Adding nodes requires a direction/tangent. As we start
        # from two images we can't do a cubic spline yet, so there are also
        # no splined tangents.
        # For the first two new nodes we use a simple tangent: the
        # unit vector pointing from the left to the right image.
        #
        # With cartesian coordinate we can use the same tangent for both
        # sides of the string.
        if self.coord_type == "cart":
            # As "get_tangent" is reimplemented in this class we call the
            # implementation of the parent ChainOfStates class.
            init_tangent = super().get_tangent(0)
            # Initial distances between left and right image
            Sk, _ = self.arc_dims
            S = Sk / (self.max_nodes+1)
            # Create first two mobile nodes
            left_step = S*init_tangent
            right_step = -S*init_tangent
        # With DLC we can't use the same tangent for both sides of the string.
        # While left_img and right_img got the same set of primitive internals,
        # they don't share the same active set U. As the DLC tangent is given
        # in the respective active set U, we have to calculate a different tanget
        # for each image.
        # TODO: Use a primitive tangent, instead of a DLC tangent.
        elif self.coord_type == "dlc":
            left_right_tangent = left_img - right_img
            l_norm = np.linalg.norm(left_right_tangent)
            Sl = l_norm / (self.max_nodes+1)
            left_step = Sl*left_right_tangent/l_norm

            right_left_tangent = right_img - left_img
            r_norm = np.linalg.norm(right_left_tangent)
            Sr = r_norm / (self.max_nodes+1)
            right_step = Sr*right_left_tangent/r_norm
        else:
            raise Exception("Invalid coord_type.")

        left_frontier = self.get_new_image(left_step, 1, 0)
        self.left_string.append(left_frontier)
        right_frontier = self.get_new_image(right_step, 2, 2)
        self.right_string.append(right_frontier)
        # cart_param_density = self.get_cur_param_density("cart")
        # coord_param_density = self.get_cur_param_density("coords")

        # Now we have four images and can calculate an initial set of tangents
        # as first derivative of the cubic spline.
        self.set_tangents()
        # The desired spacing of the nodes in the final string on the
        # normalized arclength.
        self.sk = 1 / (self.max_nodes+1)

    def get_cur_param_density(self, kind="cart"):
        if kind == "cart":
            coords = np.array([image.cart_coords for image in self.images])
            coords_ = coords.reshape(len(self.images), -1)
            diffs = coords_ - coords_[0]
        elif kind == "coords":
            image0 = self.images[0]
            diffs = np.array([image-image0 for image in self.images])
        else:
            raise Exception("Invalid kind")

        norms = np.linalg.norm(diffs, axis=1)
        # Assert that the last (rightmost) image is also the one that is the
        # farthest away from the first (leftmost) image.
        assert norms[-1] == norms.max()
        cur_param_density = norms / norms.max()
        return cur_param_density

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
    def fully_grown(self):
        """Returns wether the string is fully grown."""
        # return not ((self.left_size-1 + self.right_size-1) < self.max_nodes)
        return not (self.string_size < self.max_nodes)

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

        Tangent-calculation by splining requires the information of all
        images at once. To avoid the repeated splining of all images whenever
        a tangent is requested this method calculates all tangents and stores
        them in the self._tangents, that can be accessed via the self.tangents
        property.

        !!! Right now one must not forget to call this method
        after coordinate modification, e.g. after
        reparametrization!  Otherwise wrong (old) tangets are used. !!!
        """

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
        # Use splined tangents with cartesian coordinates
        if self.coord_type == "cart":
            return self._tangents[i]

        # With DLC we can use conventional tangents that aren't splined.

        # Upwinding tangent when the string is fully grown.
        if self.fully_grown:
            return super().get_tangent(i, kind="upwinding")

        # By definition the tangents shall point inwards during the
        # growth phase.
        cur_image = self.images[i]
        # next_ind = (i + 1) if (i <= self.lf_ind) else (i - 1)
        if i <= self.lf_ind:
            next_ind = i + 1
        else:
            next_ind = i - 1
        next_image = self.images[next_ind]
        tangent = next_image - cur_image
        tangent /= np.linalg.norm(tangent)
        return tangent

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

        np.save("gs_perp.npy", self._forces)
        return self._forces

    def reparametrize(self):
        # Non-fully-grown strings are reparametrized every cycle
        can_reparametrize = True
        self.reparam_in -= 1
        # Fully-grown strings are reparametrized only every n-th cycle
        if self.fully_grown and not (self.reparam_in == 0):
            self.log("Skipping reparametrization. Next reparametrization in "
                     f"{self.reparam_in} cycles.")
            self.set_tangents()
            return False

        # Spline displaced coordinates. 'tcks' contains all relevant information.
        # These splines will be used to interpolate all present
        # (already existing and new) images.
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
        if (not self.fully_grown) and converged(self.lf_ind):
            # Insert at the end of the left string, just before the
            # right frontier node.
            new_left_frontier = self.get_new_image(self.zero_step,
                                                   self.rf_ind, self.lf_ind)
            self.left_string.append(new_left_frontier)
            self.log("Added new left frontier node.")
        if (not self.fully_grown) and converged(self.rf_ind):
            # Insert at the end of the right string, just before the
            # current right frontier node.
            new_right_frontier = self.get_new_image(self.zero_step, self.rf_ind,
                                                    self.rf_ind)
            self.right_string.append(new_right_frontier)
            self.log("Added new right frontier node.")

        self.log(f"Current string size: {self.left_size}+{self.right_size}")
        # Reparametrize nodes
        left_inds = np.arange(self.left_size)
        right_inds = np.arange(self.max_nodes+2)[-self.right_size:]
        param_inds = np.concatenate((left_inds, right_inds))
        param_density = self.sk*param_inds
        self.log(f"New param density: " + np.array2string(param_density, precision=2))

        if self.coord_type == "cart":
            self.reparam(tcks, param_density)
            self.set_tangents()
        elif self.coord_type == "dlc":
            # coord_diffs = np.diff([image.coords for image in self.images], axis=0)
            coords_ = self.coords.reshape(len(self.images), -1)
            diffs = coords_ - coords_[0]
            norms = np.linalg.norm(diffs, axis=1)
            # Assert that the last images is also the one that is the farthest
            assert norms[-1] == norms.max()
            cur_param_density = norms / norms.max()
        else:
            raise Execption()
        self.reparam_in = self.reparam_every

        return True

    def get_additional_print(self):
        size_str = f"{self.left_size}+{self.right_size}"
        if self.fully_grown:
            size_str = "Full"
        size_info = f"String={size_str}"
        energies = np.array(self.all_energies[-1])
        barrier = (energies.max() - energies[0]) * AU2KJPERMOL
        barrier_info = f"(E_max-E_0)={barrier:.1f} kJ/mol"
        hei_ind = energies.argmax()
        hei_str = f"HEI={hei_ind+1:02d}/{energies.size:02d}"

        tot = f"Grads={self.get_image_calc_counter_sum()}"

        strs = (
            size_info,
            hei_str,
            barrier_info,
        )
        return "\t" + " ".join(strs)
