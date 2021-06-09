# [1] https://doi.org/10.1002/jcc.26495
#     Habershon, 2021

import itertools as it

import numpy as np

from pysisyphus.helpers_pure import get_molecular_radius


class HardSphereCalculator:
    def __init__(self, geom, frag_lists, kappa=1.0):
        """Intra-Image Inter-Molecular Hard-Sphere force.

        See A.2. S_2 in [1], Eqs. (A1).
        """
        self.frag_lists = frag_lists
        self.kappa = kappa

        self.frag_num = len(self.frag_lists)
        self.frag_sizes = np.array([len(frag) for frag in self.frag_lists])
        self.pair_inds = list(it.combinations(range(self.frag_num), 2))
        self.frag_inds = np.array([m for m, n in self.pair_inds])

        c3d = geom.coords3d
        frag_c3ds = [c3d[frag] for frag in self.frag_lists]
        self.frag_radii = [get_molecular_radius(frag_c3d) for frag_c3d in frag_c3ds]
        self.radii_sums = np.array(
            [self.frag_radii[i] + self.frag_radii[j] for i, j in self.pair_inds]
        )

    def get_energy(self, atoms, coords, kappa=None):
        return {"energy": 1}

    def get_forces(self, atoms, coords, kappa=None):
        if kappa is None:
            kappa = self.kappa
        c3d = coords.reshape(-1, 3)
        centroids = np.array([c3d[frag].mean(axis=0) for frag in self.frag_lists])
        mm, nn = np.triu_indices(self.frag_num, k=1)
        gdiffs = centroids[mm] - centroids[nn]
        gnorms = np.linalg.norm(gdiffs, axis=1)
        H = (gnorms < self.radii_sums).astype(int)
        N = np.zeros_like(H)
        for frag_ind, H_ in zip(self.frag_inds, H):
            N[frag_ind] += H
        nonzero = H > 0
        N *= 3 * self.frag_sizes[self.frag_inds]
        phi = np.where(N > 0, kappa / N * (gnorms - self.radii_sums), 0)
        frag_forces = phi * H * gdiffs / gnorms
        forces = np.zeros_like(c3d)
        # Distribute forces onto fragments
        for frag_ind, ff in zip(self.frag_inds, frag_forces):
            forces[self.frag_lists[frag_ind]] += ff
        # As far as I can tell so far this is actually more like a gradient.
        # as if we step along the "forces" direction the force increases
        # and ^
        forces = -forces.flatten()
        return {"energy": 1, "forces": forces}
