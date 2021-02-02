import numpy as np
from scipy.interpolate import splprep, splev

from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.intcoords.exceptions import DifferentCoordLengthsException, DifferentPrimitivesException
from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.cos.GrowingChainOfStates import GrowingChainOfStates


# [1] https://aip.scitation.org/doi/abs/10.1063/1.1691018
#     Peters, 2004
# [2] https://aip.scitation.org/doi/abs/10.1063/1.4804162
#     Zimmerman, 2013


class GrowingString(GrowingChainOfStates):

    def __init__(self, images, calc_getter, perp_thresh=0.05, param="equi",
                 reparam_every=2, reparam_every_full=3, reparam_tol=None,
                 reparam_check="rms", max_micro_cycles=5, reset_dlc=True,
                 climb=False, **kwargs):
        assert len(images) >= 2, "Need at least 2 images for GrowingString."
        if len(images) > 2:
            images = [images[0], images[-1]]
            print("More than 2 images given. Will only use first and last image!")
        if climb:
            climb = "one"
        super().__init__(images, calc_getter, climb=climb, **kwargs)

        self.perp_thresh = perp_thresh
        self.param = param
        self.reparam_every = int(reparam_every)
        self.reparam_every_full = int(reparam_every_full)
        assert self.reparam_every >= 1 and self.reparam_every_full >= 1, \
            "reparam_every and reparam_every_full must be positive integers!"
        if reparam_tol is not None:
            self.reparam_tol = float(reparam_tol)
            assert self.reparam_tol > 0
        else:
            self.reparam_tol = 1 / (self.max_nodes + 2) / 2
        self.log(f"Using reparametrization tolerance of {self.reparam_tol:.4f}")
        self.reparam_check = reparam_check
        assert self.reparam_check in ("norm", "rms")
        self.max_micro_cycles= int(max_micro_cycles)
        self.reset_dlc = bool(reset_dlc)

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
        self.new_image_inds = list()

    def get_cur_param_density(self, kind=None):
        diffs = [image - self.images[max(i-1, 0)]
                 for i, image in enumerate(self.images)]
        norms = np.linalg.norm(diffs, axis=1)
        param_density = np.cumsum(norms)
        self.log(f"Current string length={param_density[-1]:.6f}")

        # Energy weighted parametrization density
        if kind == "energy":
            prev_energies = np.array(self.all_energies[-1])

            if len(prev_energies) != len(self.images):
                return None

            mean_energies = (prev_energies[1:] + prev_energies[:-1]) / 2
            weights = mean_energies - prev_energies.min()
            # This damps everything a bit.
            weights = np.sqrt(weights)
            param_density = [0, ]
            for weight, diff in zip(weights, norms[1:]):
                assert weight > 0.
                param_density.append(param_density[-1] + weight*diff)

        param_density = np.array(param_density)
        param_density /= param_density[-1]

        return param_density

    def reset_geometries(self, ref_geometry):
        ref_typed_prims = ref_geometry.internal.typed_prims
        self.log(f"Resetting image primitives. Got {len(ref_typed_prims)} typed primitives.")
        for i in range(3):
            self.log(f"\tMicro cycle {i:d}")
            intersect = set(self.images[0].internal.typed_prims)
            for j, image in enumerate(self.images):
                image.reset_coords(ref_typed_prims)
                new_typed_prims = set(image.internal.typed_prims)
                self.log(f"\tImage {j:02d} now has {len(new_typed_prims)} typed primitives.")
                intersect = intersect & new_typed_prims

            if intersect == set(ref_typed_prims):
                ref_geometry.reset_coords(intersect)
                break
            ref_typed_prims = list(intersect)
        else:
            raise Exception("Too many reset cycles!")

    def get_new_image(self, ref_index):
        """Get new image by taking a step from self.images[ref_index] towards
        the center of the string."""
        new_img = self.images[ref_index].copy(coord_kwargs={"check_bends": True,})

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
        try:
            distance = -(new_img - tangent_img)
        except (DifferentCoordLengthsException, DifferentPrimitivesException):
            self.reset_geometries(new_img)
            distance = -(new_img - tangent_img)

        # The desired step(_length) for the new image be can be easily determined
        # from a simple rule of proportion by relating the actual distance between
        # two images to their parametrization density difference on the normalized
        # arclength and the desired spacing given by self.sk.
        #
        # Δparam_density / distance = self.sk / step
        # step = self.sk / Δparam_density * distance
        cpd = self.get_cur_param_density()
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
            chkfiles = ref_calc.get_chkfiles()
            new_img.calculator.set_chkfiles(chkfiles)
            self.log( "Set checkfiles from calculator of node "
                     f"{ref_index:02d} on calculator of new node."
            )
        except AttributeError:
            self.log("Calculator doesn't support 'get/set_chkfiles()'")
        self.images.insert(insert_ind, new_img)
        self.log(f"Created new image; inserted it before index {insert_ind}.")
        return new_img

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
    def nodes_missing(self):
        """Returns the number of nodes to be grown."""
        return (self.max_nodes + 2) - self.string_size

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

    @property
    def image_inds(self):
        return self.full_string_image_inds

    def spline(self, tangents=False):
        if (not tangents) and (self.param == "energy") and self.fully_grown:
            u = self.get_cur_param_density(kind="energy")
        else:
            u = self.get_cur_param_density()
        reshaped = self.coords.reshape(-1, self.coords_length)
        # To use splprep we have to transpose the coords.
        transp_coords = reshaped.transpose()
        # Spline in batches as scipy can't handle > 11 rows at once
        tcks, us = zip(*[splprep(transp_coords[i:i+9], s=0, k=3, u=u)
                         for i in range(0, len(transp_coords), 9)]
        )
        return tcks, us

    def reparam_cart(self, desired_param_density):
        tcks, us = self.spline()
        # Reparametrize mesh
        new_points = np.vstack([splev(desired_param_density, tck) for tck in tcks])
        # Flatten along first dimension.
        new_points = new_points.reshape(-1, len(self.images)).T
        # With a climbing image we ignore the just splined coordinates for the CI
        # and restore its original coordinates.
        for index in self.get_climbing_indices():
            new_points[index] = self.images[index].coords
            self.log(f"Skipped reparametrization of climbing image with index {index}")
        self.coords = new_points.flatten()
        # In contrast to self.reparam_dlc() we don't check if the reparametrization
        # succeeded because it can't fail ;)

    def reparam_dlc(self, desired_param_density, thresh=1e-3):
        climbing_indices = self.get_climbing_indices()
        # Reparametrization will take place along the tangent between two
        # images. The index of the tangent image depends on wether the image
        # is above or below the desired param_density on the normalized arc.
        #
        # The reparametrization is done in micro cycles, until it is converged.
        cur_param_density = self.get_cur_param_density()
        self.log(f"Density before reparametrization: {cur_param_density}")
        for i, reparam_image in enumerate(self.images[1:-1], 1):
            if i in climbing_indices:
                self.log(f"Skipped reparametrization of climbing image with index {i}")
                continue
            self.log(f"Reparametrizing node {i}")
            for j in range(self.max_micro_cycles):
                diff = (desired_param_density - cur_param_density)[i]
                self.log(f"\t{j}: Δ={diff: .6f}")
                # Do at least one pass
                if (j > 0) and (abs(diff) < thresh):
                    break
                # Negative sign: image is too far right and has to be shifted left.
                # Positive sign: image is too far left and has to be shifted right.
                sign = int(np.sign(diff))
                # Index of the tangent image. reparam_image will be shifted along
                # this direction to achieve the desired parametirzation density.
                tangent_ind = i + sign
                tangent_image = self.images[tangent_ind]
                rl = "right" if sign > 0 else "left"
                self.log(f"\t... shifting {rl} towards image {tangent_ind}")
                distance = -(reparam_image - tangent_image)

                param_dens_diff = abs(cur_param_density[tangent_ind] - cur_param_density[i])
                step_length = abs(diff) / param_dens_diff
                step = step_length * distance
                reparam_coords = reparam_image.coords + step
                reparam_image.coords = reparam_coords
                cur_param_density = self.get_cur_param_density()
            else:
                self.log(f"Reparametrization of node {i} did not converge after "
                         f"{self.max_micro_cycles} cycles. Breaking!")
                break

        cpd_str = np.array2string(cur_param_density, precision=4)
        self.log(f"Param density after reparametrization: {cpd_str}")

        # This check is disabled at it is not really applicable. While we reparametrize
        # the images the string size may vary wildly, at least in the beginning. Lets
        # say after reparametrization the distance vector between image 0 and 1 is of
        # magnitude 1 and the overall string length is 10. Then image 1 is at 0.1 w.r.t.
        # the parametrization density. If we reparametrize the remaining images the over-
        # all string size may be 8, and now image 1 suddenly sits at 1/8 = 0.125, which
        # may be already above the allowed threshold.
        # Over time the string size will equilibrate and the desired parametrization
        # density will actually be realized.
        # try:
            # # Dont check climbing images
            # np.testing.assert_allclose(
                # np.delete(cur_param_density, climbing_indices),
                # np.delete(desired_param_density, climbing_indices),
                # atol=self.reparam_tol
            # )
        # except AssertionError as err:
            # trj_str = self.as_xyz()
            # fn = "failed_reparametrization.trj"
            # with open(fn, "w") as handle:
                # handle.write(trj_str)
            # print(f"Wrote coordinates of failed reparametrization to '{fn}'")
            # raise err

        # Regenerate active set after reparametrization
        if self.reset_dlc and not self.fully_grown:
            [image.internal.set_active_set() for image in self.moving_images]
            self.log(f"Created new DLCs for {len(self.images)} string images.")
        elif self.reset_dlc:
            self.log("Skipping creation of new DLCs, as string is already fully grown.")

    def get_tangent(self, i):
        # Simple tangent, pointing at each other, for the frontier images.
        if not self.fully_grown and i in (self.lf_ind, self.rf_ind):
            next_ind = i+1 if (i <= self.lf_ind) else i-1
            tangent = self.images[next_ind] - self.images[i]
            tangent /= np.linalg.norm(tangent)
        else:
            tangent = super().get_tangent(i, kind="upwinding")

        # Converge exact mode at climbing image if requested. Use the upwinding
        # tangent as guess.
        if self.started_climbing_lanczos and (i in self.get_climbing_indices()):
            # tangent = super().get_tangent(i, kind="lanczos", lanczos_guess=tangent)
            # By using guess=None the previous Lanczos will automatically be
            # used as guess.
            tangent = super().get_tangent(i, kind="lanczos", lanczos_guess=None)
        return tangent

    @ChainOfStates.forces.getter
    def forces(self):
        if self._forces is None:
            self.calculate_forces()
        indices = range(len(self.images))
        # In constrast to NEB calculations we only use the perpendicular component
        # of the force, without any spring forces. A desired image distribution is
        # achieved via periodic reparametrization.
        perp_forces = np.array([self.get_perpendicular_forces(i) for i in indices])
        self.perp_forces_list.append(perp_forces.copy().flatten())
        # Add climbing forces
        total_forces = self.set_climbing_forces(perp_forces)
        self._forces = total_forces.flatten()
        return self._forces

    def reparametrize(self):
        reparametrized = False
        # If this counter reaches 0 reparametrization will occur.
        self.reparam_in -= 1

        self.new_image_inds = list()
        # Check if new images can be added for incomplete strings.
        if not self.fully_grown:
            perp_forces  = self.perp_forces_list[-1].reshape(len(self.images), -1)
            # Calculate norm and rms of the perpendicular force for every
            # node/image on the string.
            to_check = {
                "norm": np.linalg.norm(perp_forces, axis=1),
                "rms": np.sqrt(np.mean(perp_forces**2, axis=1)),
            }
            self.log(f"Checking frontier node convergence, threshold={self.perp_thresh:.6f}")
            # We can add a new node if the norm/rms of the perpendicular force is below
            # the threshold.
            def converged(i):
                cur_val = to_check[self.reparam_check][i]
                is_converged = cur_val <= self.perp_thresh
                conv_str = ", converged" if is_converged else ""
                self.log(f"\tnode {i:02d}: {self.reparam_check}(perp_forces)={cur_val:.6f}"
                         f"{conv_str}")
                return is_converged

            # New images are added with the same coordinates as the frontier image.
            # We force reparametrization by setting self.reparam_in to 0 to get sane
            # coordinates for the new image(s).
            if converged(self.lf_ind):
                # Insert at the end of the left string, just before the
                # right frontier node.
                new_left_frontier = self.get_new_image(self.lf_ind)
                self.new_image_inds.append(self.left_size)
                self.left_string.append(new_left_frontier)
                self.log("Added new left frontier node.")
                self.reparam_in = 0
            # If an image was just grown in the left substring the string may now
            # be fully grown, so we reavluate 'self.fully_grown' here.
            if (not self.fully_grown) and converged(self.rf_ind):
                # Insert at the end of the right string, just before the
                # current right frontier node.
                new_right_frontier = self.get_new_image(self.rf_ind)
                self.new_image_inds.append(self.left_size)
                self.right_string.append(new_right_frontier)
                self.log("Added new right frontier node.")
                self.reparam_in = 0
            self.log(f"New image indices: {self.new_image_inds}")

        self.log(
            f"Current string size is {self.left_size}+{self.right_size}="
            f"{self.string_size}. There are still {self.nodes_missing} "
            "nodes to be grown."
            if not self.fully_grown else "String is fully grown."
        )

        if self.reparam_in > 0:
            self.log("Skipping reparametrization. Next reparametrization in "
                     f"{self.reparam_in} cycles.")
        else:
            # Prepare image reparametrization
            desired_param_density = self.sk*self.full_string_image_inds
            pd_str = np.array2string(desired_param_density, precision=4)
            self.log(f"Desired param density: {pd_str}")

            # Reparametrize images.
            if self.coord_type == "cart":
                self.reparam_cart(desired_param_density)
            elif self.coord_type == "dlc":
                self.reparam_dlc(desired_param_density, thresh=self.reparam_tol)
            else:
                raise Exception("How did you get here?")

            self.reparam_in = self.reparam_every_full if self.fully_grown \
                              else self.reparam_every
            reparametrized = True
            with open("reparametrized.trj", "w") as handle:
                handle.write(self.as_xyz())

        return reparametrized

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

        strs = (
            size_info,
            hei_str,
            barrier_info,
        )
        return "\t" + " ".join(strs)
