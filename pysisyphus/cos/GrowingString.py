#/!usr/bin/env python3

# See [1] 10.1063/1.1691018

from copy import copy

import numpy as np
from scipy.interpolate import splprep, splev

from pysisyphus.cos.GrowingChainOfStates import GrowingChainOfStates


class GrowingString(GrowingChainOfStates):

    def __init__(self, images, calc_getter, **kwargs):
        assert len(images) == 2, "Can only start from 2 images for now"
        super().__init__(images, calc_getter, **kwargs)

        left_img, right_img = self.images

        self.left_string = [left_img, ]
        self.right_string = [right_img, ]

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
        add_ks = np.arange(self.max_nodes)

        self.points = [self.images[0].coords]
        self.conv_points = [self.images[0].coords]

        # To add nodes we need a direction/tangent. As we start
        # from two images we can't spline yet, so we got no splined
        # tangents.
        # Instead we use the unit vector pointing from the
        # left image to the right image.
        init_tangent = super().get_tangent(0)
        # S = .75
        Sk, _ = self.arc_dims
        S = Sk / (self.max_nodes+1)
        print("S", S)
        # Create first two mobile nodes
        left_img, right_img = self.images
        # new_left_coords = left_img.coords + sk*init_tangent
        # new_right_coords = right_img.coords - sk*init_tangent
        new_left_coords = left_img.coords + S*init_tangent
        new_right_coords = right_img.coords - S*init_tangent
        left_frontier = self.get_new_image(new_left_coords, 1)
        self.left_string.append(left_frontier)
        right_frontier = self.get_new_image(new_right_coords, 2)
        self.right_string.append(right_frontier)

        def ind_func(perp_force, tol=0.5):
            return int(np.linalg.norm(perp_force) <= tol)

        tcks, us = self.spline()
        print("us", us)
        for i in range(35):
            Sk, cur_mesh = self.arc_dims
            print("cur_mesh", cur_mesh)
            a1, a2 = cur_mesh[self.lf_ind], cur_mesh[self.lf_ind+1]
            # Step length on the normalized arclength
            S = Sk / (self.max_nodes+1)
            sk = S / Sk
            print(f"cycle {i}, total arclength Sk={Sk:.4f}, sk={sk:.4f}")
            forces = self.forces
            # print("forces")
            # print(forces.reshape(-1, 3))
            # print()
            self.cur_forces = forces.reshape(-1, 3)
            tangents = np.vstack([splev(cur_mesh, tck, der=1) for tck in tcks]).T
            norms = np.linalg.norm(tangents, axis=1)
            tangents = tangents / norms[:,None]
            right = len(self.left_string)
            # Right tangents shall point inward
            tangents[right:] *= -1
            self.tangents = tangents
            # print("tangents")
            # print(tangents)
            # print()
            # Perp forces = org. forces - parallel part
            tf = tangents.flatten()
            perp_forces = list()
            force_per_img = forces.reshape(-1, 3)
            for i, (force, tang) in enumerate(zip(force_per_img, tangents)):
                pf = force - force.dot(tang)*tang
                # print(i, pf)
                perp_forces.append(pf)
            perp_forces = np.array(perp_forces).flatten()
            self.perp_forces = perp_forces.reshape(-1, 3)
            perp_norms = np.linalg.norm(self.perp_forces, axis=1)
            print("perp_norms")
            print(np.array2string(perp_norms, precision=2))
            step = 0.05 * perp_forces
            step_norms = np.linalg.norm(step.reshape(-1, 3), axis=1)
            # print("step_norms")
            # print(step_norms)
            # print()

            self.coords += step
            tcks, us = self.spline()

            print(f"a1={a1}, a2={a2}")
            # a1 += ind_func(self.perp_forces[self.lf_ind])*sk
            # a2 -= ind_func(self.perp_forces[self.rf_ind])*sk
            if self.images_left and ind_func(self.perp_forces[self.lf_ind]):
                # Insert at the end of the left string, just before the 
                # right frontier node.
                new_left_frontier = self.get_new_image(self.dummy_coords, self.rf_ind)
                self.left_string.append(new_left_frontier)
                a1 += sk
                # print("perp force norm on left frontier is below threshold!")
                # print("propagate left, a1")
                print("adding new left frontier node.")
            if self.images_left and ind_func(self.perp_forces[self.rf_ind]):
                # Insert at the end of the right string, just before the 
                # current right frontier node.
                new_right_frontier = self.get_new_image(self.dummy_coords, self.rf_ind)
                self.right_string.append(new_right_frontier)
                a2 -= sk
                print("adding new right frontier node.")
                # print("perp force norm on right frontier is below threshold!")
                # print("propagate right, a2")
            print("string size:", self.left_size, "+", self.right_size)
            # print(f"got {self.left_size} images in left string")
            # print(f"got {self.right_size} images in right string")
            # print(f"a1={a1}, a2={a2}")
            left_density = np.linspace(0, a1, self.left_size)
            right_density = np.linspace(a2, 1, self.right_size)
            param_density = np.concatenate((left_density, right_density))
            print("new param density", param_density)
            self.reparam(tcks, param_density)
            print()


        self.conv_points = [i.coords for i in self.images]
        # import pdb; pdb.set_trace()
        # print("diffs", np.diff(self.conv_points))
