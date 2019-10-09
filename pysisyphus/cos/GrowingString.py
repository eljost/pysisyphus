#/!usr/bin/env python3

# See [1] 10.1063/1.1691018

import numpy as np
from scipy.interpolate import splprep, splev

from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.cos.GrowingChainOfStates import GrowingChainOfStates


class GrowingString(GrowingChainOfStates):

    def __init__(self, images, calc_getter, perp_thresh=0.05,
                 reparam_every=3, reparam_tol=None, **kwargs):
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
        if reparam_tol is not None:
            self.reparam_tol = float(reparam_tol)
            assert self.reparam_tol > 0
        else:
            self.reparam_tol = 1 / (self.max_nodes + 2) / 2
        self.log(f"Using reparametrization tolerance of {self.reparam_tol:.4e}")

        left_img, right_img = self.images

        self.left_string = [left_img, ]
        self.right_string = [right_img, ]

        # The desired spacing of the nodes in the final string on the
        # normalized arclength.
        self.sk = 1 / (self.max_nodes+1)

        self.reparam_in = reparam_every
        self._tangents = None
        self.tangent_list = list()
        self.perp_forces_list = list()
        self.coords_list = list()

        left_frontier = self.get_new_image(self.lf_ind)
        self.left_string.append(left_frontier)
        right_frontier = self.get_new_image(self.rf_ind)
        self.right_string.append(right_frontier)

        if self.coord_type == "cart":
            self.set_tangents()

    def get_cur_param_density(self, kind="cart"):
        if kind == "cart":
            coords = np.array([image.cart_coords for image in self.images])
            coords_ = coords.reshape(len(self.images), -1)
            diffs = coords_ - coords_[0]
        elif kind == "coords":
            image0 = self.images[0]
            # This way, even with DLC all differences will be given in the
            # active set of image0.
            diffs = np.array([image0-image for image in self.images])
        else:
            raise Exception("Invalid kind")

        norms = np.linalg.norm(diffs, axis=1)
        cur_param_density = norms / norms.max()
        # Assert that the last (rightmost) image is also the one that is the
        # farthest away from the first (leftmost) image.
        assert norms[-1] == norms.max(), \
            "Unexpected parametrization density. Expected the last " \
            "(rightmost) image to be the farthest image, but this is " \
            "not the case. Current parametrization density is: " \
           f"{cur_param_density}."
        return cur_param_density

    def get_new_image(self, ref_index):
        """Get new image by taking a step from self.images[ref_index] towards
        the center of the string."""
        new_img = self.images[ref_index].copy()

        if ref_index <= self.lf_ind:
            tangent_ind = ref_index + 1
            insert_ind = tangent_ind
        else:
            tangent_ind = ref_index - 1
            insert_ind = ref_index
        tangent_img = self.images[tangent_ind]

        # (new_img - tangent_img) points from tangent_img towards new_img.
        # As we want to derive a new image from new_img, we have to step
        # against this vector, so we have to multiply by -1.
        # Why don't we just use (tangent_img - new_img) to get the right
        # direction? In DLC the resulting distance would then be given in
        # the active set U of tangent_img, but we need it in the active set U
        # of new_img.
        # Formulated the other way around the same expression can be used for
        # all coord types.
        distance = -(new_img - tangent_img)

        # The desired step(_length) for the new image be can be easily determined
        # from a simple rule of proportion by relating the actual distance between
        # two images to their parametrization density difference on the normalized
        # arclength and the desired spacing given by self.sk.
        #
        # Δparam_density / distance = self.sk / step
        # step = self.sk / Δparam_density * distance
        cpd = self.get_cur_param_density("coords")
        # As we always want to step in the direction of 'distance' we just take
        # the absolute value of the difference, as we are not interested in the
        # sign.
        param_dens_diff = abs(cpd[ref_index] - cpd[tangent_ind])
        step_length = self.sk / param_dens_diff
        step = step_length * distance

        new_coords = new_img.coords + step
        new_img.coords = new_coords
        new_img.set_calculator(self.calc_getter())
        ref_calc = self.images[ref_index].calculator
        try:
            ref_calc.propagate_wavefunction(new_img.calculator)
            self.log( "Set wavefunction data from calculator of node "
                     f"{ref_index:02d} on calculator of new node."
            )
        except AttributeError:
            self.log("Calculator doesn't support 'propagte_wavefunction()'")
        self.images.insert(insert_ind, new_img)
        self.log(f"Created new image; inserted it before index {insert_ind}.")
        return new_img

        # self.images.insert(insert_ind, new_img)
        # # Take smaller steps, as the internal-cartesian-backconversion may be
        # # unstable for bigger steps.
        # steps = 10
        # step = step_length * distance/steps
        # for i in range(steps):
            # new_coords = new_img.coords + step
            # new_img.coords = new_coords
            # cpd = self.get_cur_param_density("coords")
            # try:
                # if new_img.internal.backtransform_failed:
                    # import pdb; pdb.set_trace()
            # except AttributeError:
                # pass
            # print(f"{i:02d}: {cpd}")

        # # self.images.insert(insert_ind, new_img)
        # # self.log(f"Created new image; inserted it before index {insert_ind}.")

        # cpd = self.get_cur_param_density("coords")
        # self.log(f"Current param_density: {cpd}")

        # return new_img

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
        """Returns wether the string is fully grown. Don't count the first
        and last node."""
        return not ((self.string_size - 2) < self.max_nodes)

    @property
    def lf_ind(self):
        """Index of the left frontier node in self.images."""
        return len(self.left_string)-1

    @property
    def rf_ind(self):
        """Index of the right frontier node in self.images."""
        return self.lf_ind+1

    @property
    def full_string_image_inds(self):
        left_inds = np.arange(self.left_size)
        right_inds = np.arange(self.max_nodes+2)[-self.right_size:]
        image_inds = np.concatenate((left_inds, right_inds))
        return image_inds

    def spline(self):
        reshaped = self.coords.reshape(-1, self.coords_length)
        # To use splprep we have to transpose the coords.
        transp_coords = reshaped.transpose()
        # Spline in batches as scipy can't handle > 11 rows at once
        tcks, us = zip(*[splprep(transp_coords[i:i+9], s=0, k=3)
                         for i in range(0, len(transp_coords), 9)]
        )
        return tcks, us

    def reparam_cart(self, tcks, param_density):
        # Reparametrize mesh
        new_points = np.vstack([splev(param_density, tck) for tck in tcks])
        # Flatten along first dimension.
        new_points = new_points.reshape(-1, len(self.images))
        self.coords = new_points.transpose().flatten()

    # def reparam_dlc(self, cur_param_density, desired_param_density, thresh=1e-3):
        # # Reparametrization will take place along the tangent between two
        # # images. The index of the tangent image depends on wether the image
        # # is above or below the desired param_density on the normalized arc.
        # diffs = desired_param_density - cur_param_density
        # # Negative sign: image is too far right and has to be shifted left.
        # # Positive sign: image is too far left and has to be shifted right.
        # signs = np.sign(diffs).astype(int)
        # # TODO: multiple passes of this loop to get a tighter convergence,
        # # so a lower atol can used in the np.testing method.
        # for i, (diff, sign) in enumerate(zip(diffs, signs)):
            # if abs(diff) < thresh:
                # continue
            # reparam_image = self.images[i]
            # # Index of the tangent image. reparam_image will be shifted along
            # # this direction to achieve the desired parametirzation density.
            # tangent_ind = i + sign
            # tangent_image = self.images[tangent_ind]
            # distance = -(reparam_image - tangent_image)

            # param_dens_diff = abs(cur_param_density[tangent_ind] - cur_param_density[i])
            # step_length = abs(diff) / param_dens_diff
            # step = step_length * distance
            # reparam_coords = reparam_image.coords + step
            # reparam_image.coords = reparam_coords
            # cur_param_density = self.get_cur_param_density("coords")
        # np.testing.assert_allclose(cur_param_density, desired_param_density,
                                   # # atol=max(5e-2, 5*thresh))
                                   # atol=thresh)

        # # Regenerate active set after reparametrization
        # # [image.internal.set_active_set() for image in self.moving_images]

    def reparam_dlc(self, cur_param_density, desired_param_density, thresh=1e-3):
        # Reparametrization will take place along the tangent between two
        # images. The index of the tangent image depends on wether the image
        # is above or below the desired param_density on the normalized arc.

        # This implementation assumes that the reparametrization step take is not
        # too big, so the internal-cartesian-transformation doesn't fail.
        # Adding new images is done with smaller steps to avoid this problem.
        # As every images is added only once, but may be reparametrized quite often
        # we try to do the reparametrization in one step.
        # A safer approach would be to do it in multiple smaller steps.

        for i, reparam_image in enumerate(self.images[1:-1], 1):
            diff = (desired_param_density - cur_param_density)[i]
            # Negative sign: image is too far right and has to be shifted left.
            # Positive sign: image is too far left and has to be shifted right.
            sign = int(np.sign(diff))
            if abs(diff) < thresh:
                continue
            # Index of the tangent image. reparam_image will be shifted along
            # this direction to achieve the desired parametirzation density.
            tangent_ind = i + sign
            tangent_image = self.images[tangent_ind]
            distance = -(reparam_image - tangent_image)

            param_dens_diff = abs(cur_param_density[tangent_ind] - cur_param_density[i])
            step_length = abs(diff) / param_dens_diff
            step = step_length * distance
            reparam_coords = reparam_image.coords + step
            reparam_image.coords = reparam_coords
            cur_param_density = self.get_cur_param_density("coords")
        self.log(f"Current param density: {cur_param_density}")
        np.testing.assert_allclose(cur_param_density, desired_param_density,
                                   atol=self.reparam_tol)

        # Regenerate active set after reparametrization
        # [image.internal.set_active_set() for image in self.moving_images]

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

        # With DLC we can use conventional tangents that can be calculated
        # without splining.

        # Upwinding tangent when the string is fully grown.
        if self.fully_grown:
            return super().get_tangent(i, kind="upwinding")

        # During the growth phase we use simple tangents that always point
        # towards the center of the string.
        cur_image = self.images[i]
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

        # We can add new nodes if the string is not yet fully grown
        # and if the frontier nodes are converged below self.perp_thresh.
        # Right now we add the new image(s) with a zero step, so they got
        # the same coordinates as the respective frontier geometry.
        # We then rely on reparametrization to assign the correct coordinates.
        if (not self.fully_grown) and converged(self.lf_ind):
            # Insert at the end of the left string, just before the
            # right frontier node.
            new_left_frontier = self.get_new_image(self.lf_ind)
            self.left_string.append(new_left_frontier)
            self.log("Added new left frontier node.")
        if (not self.fully_grown) and converged(self.rf_ind):
            # Insert at the end of the right string, just before the
            # current right frontier node.
            new_right_frontier = self.get_new_image(self.rf_ind)
            self.right_string.append(new_right_frontier)
            self.log("Added new right frontier node.")

        self.log(f"Current string size: {self.left_size}+{self.right_size}")

        # Prepare node reparametrization
        desired_param_density = self.sk*self.full_string_image_inds
        pd_str = np.array2string(desired_param_density, precision=3)
        self.log(f"Desired param density: {pd_str}")

        # TODO: Add some kind of threshold and only reparametrize when
        # the deviation from the desired param_density is above the threshold.
        if self.coord_type == "cart":
            self.reparam_cart(tcks, desired_param_density)
            self.set_tangents()
        elif self.coord_type == "dlc":
            cur_param_density = self.get_cur_param_density("coords")
            self.reparam_dlc(cur_param_density, desired_param_density)
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
