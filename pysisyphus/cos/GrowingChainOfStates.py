#!/usr/bin/env python3
from copy import copy

import numpy as np
from scipy.interpolate import splprep, splev

from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.Geometry import Geometry


class GrowingChainOfStates(ChainOfStates):

    def __init__(self, images, calc_getter, max_nodes=10,
                 **kwargs):
        super().__init__(images, **kwargs)

        self.max_nodes = max_nodes
        self.calc_getter = calc_getter
        self.atoms = copy(self.images[0].atoms)
    
    def get_new_image(self, coords, index):
        new_image = Geometry(self.atoms, coords)
        new_image.set_calculator(self.calc_getter())
        self.images.insert(index, new_image)
        self.log(f"Create new image; insert it before index {index}.")
        return new_image

    @property
    def dummy_coords(self):
        return np.empty_like(self.images[0].coords)

    @property
    def arc_dims(self):
        coords = [image.coords for image in self.images]
        cds = [0, ]
        for i in range(len(coords)-1):
            diff = np.linalg.norm(coords[i+1]-coords[i])
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
        new_node = Geometry(self.atoms, new_coords)
        new_node.set_calculator(self.calc_getter())
        self.images.insert(k+1, new_node)
        return new_node
