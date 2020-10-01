#/!usr/bin/env python3

# See [1] 10.1063/1.1691018

import numpy as np
from scipy.interpolate import splprep, splev

from pysisyphus.cos.GrowingChainOfStates import GrowingChainOfStates


class GrowingString(GrowingChainOfStates):

    def __init__(self, images, calc_getter, max_cycles=75, **kwargs):
        assert len(images) >= 2, "Need at least 2 images for GrowingString."
        if len(images) >= 2:
            images = [images[0], images[-1]]
            print("More than 2 images were supplied! Will only use the "
                  "first and last images to start the GrowingString."
            )
        super().__init__(images, calc_getter, **kwargs)

        self.max_cycles = max_cycles
        left_img, right_img = self.images

        self.left_string = [left_img, ]
        self.right_string = [right_img, ]

        self.tangent_list = list()
        self.perp_force_list = list()
        self.coords_list = list()

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

    def run(self):
        # To add nodes we need a direction/tangent. As we start
        # from two images we can't do a cubic spline yet, so we
        # can't use splined tangents.
        # To start off we use a "classic" tangent instead, that
        # is the unit vector pointing from the left image, to
        # the right image.
        init_tangent = super().get_tangent(0)
        Sk, _ = self.arc_dims
        S = Sk / (self.max_nodes+1)
        # Create first two mobile nodes
        new_left_coords = S*init_tangent
        new_right_coords = - S*init_tangent
        left_frontier = self.get_new_image(new_left_coords, 1, 0)
        self.left_string.append(left_frontier)
        right_frontier = self.get_new_image(new_right_coords, 2, 2)
        self.right_string.append(right_frontier)

        def ind_func(perp_force, tol=0.5):  # lgtm [py/unused-local-variable]
            return int(np.linalg.norm(perp_force) <= tol)

        # Step length on the normalized arclength
        sk = 1 / (self.max_nodes+1)  # lgtm [py/unused-local-variable]
        tcks, us = self.spline()
        for self.cur_cycle in range(self.max_cycles):
            Sk, cur_mesh = self.arc_dims
            print(f"Cycle {self.cur_cycle:03d}, total arclength Sk={Sk:.4f}")
            forces = self.forces
            self.cur_forces = forces.reshape(-1, 3)
            tangents = np.vstack([splev(cur_mesh, tck, der=1) for tck in tcks]).T
            norms = np.linalg.norm(tangents, axis=1)
            tangents = tangents / norms[:,None]
            # Tangents of the right string shall point towards the center.
            tangents[self.rf_ind:] *= -1
            self.tangents = tangents
            self.tangent_list.append(self.tangents)

            # Perpendicular force component
            # Dot product between rows in one line
            # np.einsum("ij,ij->i", tangents, forces.reshape(-1, 3))
            force_per_img = forces.reshape(-1, 3)
            self.perp_forces = np.array(
                [force - force.dot(tangent)*tangent
                 for force, tangent in zip(force_per_img, tangents)]
            )
            self.perp_force_list.append(self.perp_forces)
            np.save("gs_ref_perp.npy", self.perp_forces)

            # Take step
            step = 0.05 * self.perp_forces.flatten()
            step_norms = np.linalg.norm(step.reshape(-1, 3), axis=1)
            self.coords += step
            self.coords_list.append(self.coords)
            # Spline displaced coordinates
            tcks, us = self.spline()

            # Check if we can add new nodes
            if self.images_left and ind_func(self.perp_forces[self.lf_ind]):
                # Insert at the end of the left string, just before the 
                # right frontier node.
                new_left_frontier = self.get_new_image(self.zero_step, self.rf_ind, self.lf_ind)
                self.left_string.append(new_left_frontier)
                print("Added new left frontier node.")
            if self.images_left and ind_func(self.perp_forces[self.rf_ind]):
                # Insert at the end of the right string, just before the 
                # current right frontier node.
                new_right_frontier = self.get_new_image(self.zero_step, self.rf_ind, self.rf_ind)
                self.right_string.append(new_right_frontier)
                print("Added new right frontier node.")

            print("Current string size:", self.left_size, "+", self.right_size)
            # Reparametrize nodes
            left_inds = np.arange(self.left_size)
            right_inds = np.arange(self.max_nodes+2)[-self.right_size:]
            param_inds = np.concatenate((left_inds, right_inds))
            param_density = sk*param_inds
            print("New param density", np.array2string(param_density, precision=2))
            self.reparam(tcks, param_density)
