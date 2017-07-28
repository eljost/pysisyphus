#!/usr/bin/env python3

import numpy as np

from cos.ChainOfStates import ChainOfStates

# [1] http://aip.scitation.org/doi/pdf/10.1063/1.1323224
# 
# https://github.com/cstein/neb/blob/master/neb/neb.py


class NEB(ChainOfStates):

    def __init__(self, images):
        super(NEB, self).__init__(images)

        self.perp_forces = list()
        self.par_forces = list()

    @property
    def forces(self):
        inner_indices = list(range(1, len(self.images)-1))
        [self.images[i].calc_energy_and_forces() for i in inner_indices]
        self.fix_endpoints()
        
        self.make_tangents()

        perp_forces = np.array([self.get_perpendicular_forces(i)
                                for i in inner_indices])
        par_forces = np.array([self.get_parallel_forces(i)
                               for i in inner_indices])
        total_forces = perp_forces + par_forces

        # Store the calculated forces so they can be used later on, e.g.
        # for plotting. total_forces is saved in self._forces.
        self.perp_forces.append(perp_forces)
        self.par_forces.append(par_forces)

        for ii, tf in zip(inner_indices, total_forces):
            self.images[ii].forces = tf
        # Here we also use the zero forces of educt and product
        self._forces  = np.concatenate([image.forces for image in self.images])
        return self._forces
