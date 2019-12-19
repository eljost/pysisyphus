#!/usr/bin/env python3

import numpy as np

from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.Geometry import Geometry


class GrowingChainOfStates(ChainOfStates):

    def __init__(self, images, calc_getter, max_nodes=10,
                 **kwargs):
        super().__init__(images, **kwargs)

        self.max_nodes = max_nodes
        self.calc_getter = calc_getter
        self.zero_step = np.zeros_like(self.images[0].coords)

    # TODO: remove this as it is seems to be never used?!
    # def get_new_image(self, step, before_index, ref_index):
        # new_image = self.images[ref_index].copy()
        # new_coords = new_image.coords + step
        # new_image.coords = new_coords
        # new_image.set_calculator(self.calc_getter())
        # self.images.insert(before_index, new_image)
        # self.log(f"Created new image; inserted it before index {before_index}.")
        # return new_image

    def get_new_image_from_coords(self, coords, index):
        new_image = Geometry(self.image_atoms, coords,
                             coord_type=self.coord_type,
                             prim_indices=self.prim_indices)
        new_image.set_calculator(self.calc_getter())
        self.images.insert(index, new_image)
        self.log(f"Create new image; insert it before index {index}.")
        return new_image

    @property
    def arc_dims(self):
        cds = [0, ]
        for i, image in enumerate(self.images[:-1]):
            next_image = self.images[i+1]
            diff = np.linalg.norm(next_image - image)
            cds.append(diff)
        cds = np.cumsum(cds)
        tot_length = cds[-1]
        norm_cds = cds / cds.max()
        return tot_length, norm_cds

    def new_node_coords(self, k):
        l = (self.max_nodes-k) / (self.max_nodes+1-k)
        kth_coords = self.images[k].coords
        last_coords = self.images[-1].coords
        new_coords = l*kth_coords + (1-l)*last_coords
        return new_coords

    def set_new_node(self, k):
        new_coords = self.new_node_coords(k)
        new_node = Geometry(self.image_atoms, new_coords)
        new_node.set_calculator(self.calc_getter())
        self.images.insert(k+1, new_node)
        return new_node

    def prepare_opt_cycle(self, *args, **kwargs):
        parent_result = super().prepare_opt_cycle(*args, **kwargs)

        # Compare size of coords arrays to determine if new nodes
        # were added in the last reparametrization.
        last_size = self.coords_list[-1].size
        length_changed = last_size != self.coords.size
        return parent_result or length_changed
