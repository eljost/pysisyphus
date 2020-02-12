#!/usr/bin/env python3

import numpy as np

from pysisyphus.cos.ChainOfStates import ChainOfStates

# [1] http://aip.scitation.org/doi/pdf/10.1063/1.1323224
#     10.1063/1.1323224
# [2] http://onlinelibrary.wiley.com/doi/10.1002/jcc.20780/pdf
#     10.1002/jcc.20780
# [3] https://aip.scitation.org/doi/10.1063/1.2841941
#     Sheppard, 2008
# [4] https://aip.scitation.org/doi/pdf/10.1063/1.1636455
#     Trygubenko, 2004
# [5] Nudged Elastic Band Method for Finding Minimum Energy Paths of Transitions
#     Hannes Jónsson , Greg Mills , Karsten W. Jacobsen
# https://github.com/cstein/neb/blob/master/neb/neb.py


class NEB(ChainOfStates):

    def __init__(self, images, variable_springs=False,
                 k_max=0.3, k_min=0.1, perp_spring_forces=None, **kwargs):
        super(NEB, self).__init__(images, **kwargs)

        assert(k_max > k_min), "k_max has to be bigger than k_min!"
        self.variable_springs = variable_springs
        self.k_max = k_max
        self.k_min = k_min
        self.perp_spring_forces = perp_spring_forces

        self.delta_k = self.k_max - self.k_min
        self.k = list()

    def update_springs(self):
        # Check if there are enough springs
        if (len(self.k) != len(self.images)-1):
            self.k = np.full(len(self.images)-1, self.k_min)
        if self.variable_springs:
            self.set_variable_springs()

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

    def get_spring_forces(self, i):
        if i not in self.moving_indices:
            return self.zero_vec.copy()

        if (i == 0) or (i == len(self.images) - 1):
            # We can't use the last image index because there is one
            # spring less than there are images.
            spring_index = min(i, len(self.images)-2)
            return self.k[spring_index] * self.get_tangent(i)

        prev_coords = self.images[i-1].coords
        ith_coords = self.images[i].coords
        next_coords = self.images[i+1].coords
        spring_forces = (
            self.k[i] * (next_coords-ith_coords) - (ith_coords-prev_coords)
        )
        return spring_forces

    def get_quenched_dneb_forces(self, i):
        """See [3], Sec. VI and [4] Sec. D."""
        if not self.perp_spring_forces or (i not in self.moving_indices):
            return self.zero_vec.copy()
        forces = self.images[i].forces
        tangent = self.get_tangent(i)
        perp_forces = forces - forces.dot(tangent)*tangent
        spring_forces = self.get_spring_forces(i)
        tangent = self.get_tangent(i)
        perp_spring_forces = spring_forces - spring_forces.dot(tangent)*tangent
        dneb_forces = perp_spring_forces - perp_spring_forces.dot(perp_forces)*perp_forces
        perp_norm = np.linalg.norm(perp_forces)
        perp_spring_norm = np.linalg.norm(perp_spring_forces)

        # Switching function to quench the dneb forces
        # Eq. (15) in [3]
        #
        # If norm(perp_force) >> norm(perp_spring_forces): dneb_factor ~ 1
        # If norm(perp_force) << norm(perp_spring_forces): dneb_factor ~ 0
        #
        # If the perpendicular spring force is much bigger than the
        # perpendicular force the DNEB forces is nearly fully quenched.
        dneb_factor = 2/np.pi * np.arctan2(perp_norm**2, perp_spring_norm**2)
        dneb_forces_quenched = dneb_factor * dneb_forces

        # An alternative switchting function is given in [5], Eq. (10)
        # f(phi) = 1/2 * (1 + cos(pi*cos(theta)))
        # f -> 0 for a straight path (theta -> 0°)
        # f -> 1 for a perpendicular path (theta -> 90°)
        # cos(theta) = (R_(i+1) - R_i) * (R_i - R_(i-1)) / (norm of numerator)

        return dneb_forces_quenched

    def get_parallel_forces(self, i):
        if i not in self.moving_indices:
            return self.zero_vec.copy()

        if (i == 0) or (i == len(self.images) - 1):
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

        self.update_springs()

        indices = range(len(self.images))
        total_forces = np.array(
            [self.get_parallel_forces(i) + self.get_perpendicular_forces(i)
             + self.get_quenched_dneb_forces(i)
             for i in indices]
        )
        total_forces = self.set_climbing_forces(total_forces)
        self._forces = np.array(total_forces).flatten()
        return self._forces
