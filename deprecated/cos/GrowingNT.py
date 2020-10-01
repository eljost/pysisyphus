#/!usr/bin/env python3

# See [1] 10.1063/1.1885467

from copy import copy

import numpy as np

from pysisyphus.cos.ChainOfStates import ChainOfStates
from pysisyphus.Geometry import Geometry


class GrowingNT(ChainOfStates):

    def __init__(self, images, calc_getter, eps, damp, max_nodes=10,
                 readjust=True, **kwargs):
        super().__init__(images, **kwargs)

        self.calc_getter = calc_getter
        self.max_nodes = max_nodes
        self.eps = eps
        self.damp = damp
        self.readjust = readjust

        self.I = np.eye(self.images[0].coords.size)

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
        # print(f"made node {k+1}")
        return new_node

    def run(self):
        # add_ks = np.arange(self.max_nodes-len(self.images))
        add_ks = np.arange(self.max_nodes)
        # import pdb; pdb.set_trace()
        # Initial rStart at 
        # grad0 = self.images[0].gradient
        # norm_grad0 = grad0 / np.linalg.norm(grad0)
        r = self.get_tangent(0)[:,None]
        # r = norm_grad0[:,None]

        self.points = [self.images[0].coords]
        self.conv_points = [self.images[0].coords]
        readjusted = False
        for k in add_ks:
            new_node = self.set_new_node(k)
            if self.readjust and (not readjusted) and k > self.max_nodes / 2:
                # Adapt search direction
                print("old search direction", r)
                r = self.get_tangent(k+1)[:,None]
                print("adapted search direction", r)
                readjusted = True
            # print("new node", new_node.coords)
            Dr = r.dot(r.T)
            Pr = self.I - Dr
            # Do projections and correction
            for i in range(45):
                grad = new_node.gradient
                proj_grad = Pr.dot(grad)
                norm = np.linalg.norm(proj_grad)
                # print(f"cycle {i} norm is {norm:.3f}")
                if norm < self.eps:
                    print(f"{i:02d} microcycles, norm={norm:.4f}")
                    break
                p = -self.damp*proj_grad
                new_node.coords += p
                self.points.append(new_node.coords.copy())
            self.conv_points.append(new_node.coords.copy())
        self.conv_points.append(self.images[-1].coords.copy())
