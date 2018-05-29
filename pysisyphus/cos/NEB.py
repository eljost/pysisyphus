#!/usr/bin/env python3

import numpy as np

from pysisyphus.cos.ChainOfStates import ChainOfStates

# [1] http://aip.scitation.org/doi/pdf/10.1063/1.1323224
#     10.1063/1.1323224
# [2] http://onlinelibrary.wiley.com/doi/10.1002/jcc.20780/pdf
#     10.1002/jcc.20780
# https://github.com/cstein/neb/blob/master/neb/neb.py


class NEB(ChainOfStates):

    def __init__(self, images, variable_springs=False,
                 k_max=0.3, k_min=0.01, **kwargs):
        super(NEB, self).__init__(images, **kwargs)

        assert(k_max > k_min), "k_max has to be bigger than k_min!"
        self.variable_springs = variable_springs
        self.k_max = k_max
        self.k_min = k_min

        self.delta_k = self.k_max - self.k_min
        self.k = list()

        self.perp_forces = list()
        self.par_forces = list()

    def update_springs(self):
        self.k = np.full(len(self.images)-1, self.k_min)

    def set_variable_springs(self):
        shifted_energies = self.energy - self.energy.min()
        energy_max = max(shifted_energies)
        energy_ref = .85 * energy_max
        for i in range(len(self.k)):
            # The ith spring connects images i-1 and i.
            e_i = i + 1
            ith_energy = max(shifted_energies[e_i], shifted_energies[e_i-1])
            if ith_energy < energy_ref:
                self.k[i] = self.k_min
            else:
                self.k[i] = (self.k_max - self.delta_k
                             * (energy_max - ith_energy)
                             / (energy_max - energy_ref)
                )
        self.log("updated springs: " + self.fmt_k())

    def fmt_k(self):
        return ", ".join([str(f"{k:.03f}") for k in self.k])

    @property
    def parallel_forces(self):
        indices = range(len(self.images))
        par_forces = [self.get_parallel_forces(i) for i in indices]
        return np.array(par_forces).flatten()

    def get_parallel_forces(self, i):
        # Check if there are enough springs
        if (len(self.k) is not len(self.images)-1):
            self.update_springs()

        if i not in self.moving_indices:
            return self.zero_vec

        if (i is 0) or (i is len(self.images) - 1):
            # We can't use the last image index because there is one
            # spring less than there are images.
            spring_index = min(i, len(self.images)-2)
            return self.k[spring_index] * self.get_tangent(i)

        prev_coords = self.images[i-1].coords
        ith_coords = self.images[i].coords
        next_coords = self.images[i+1].coords
        return (self.k[i] * (np.linalg.norm(next_coords-ith_coords)
                             - np.linalg.norm(ith_coords-prev_coords)
               ) * self.get_tangent(i)
        )

    # See https://stackoverflow.com/a/15786149
    # This way we can reuse the parents setter.
    @ChainOfStates.forces.getter
    def forces(self):
        if self._forces is None:
            self.calculate_forces()

        indices = range(len(self.images))
        total_forces = np.array(
            [self.get_parallel_forces(i) + self.get_perpendicular_forces(i)
             for i in indices]
        )
        total_forces = self.set_climbing_forces(total_forces)
        self._forces = np.array(total_forces).flatten()
        return self._forces
