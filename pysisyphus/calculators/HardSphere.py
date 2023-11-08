# [1] https://doi.org/10.1002/jcc.26495
#     Habershon, 2021

import itertools as it

import numpy as np

from pysisyphus.helpers_pure import get_molecular_radius


class HardSphere:
    def __init__(
        self,
        geom,
        frags,
        kappa=1.0,
        permutations=False,
        frag_radii=None,
        radii_offset=0.9452,
    ):
        """Intra-Image Inter-Molecular Hard-Sphere force.

        See A.2. in [1], Eq. (A1).
        """
        self.frags = frags
        self.kappa = kappa

        self.frag_num = len(self.frags)
        self.frag_sizes = np.array([len(frag) for frag in self.frags])
        self.pair_inds = np.array(list(it.combinations(range(self.frag_num), 2)))
        it_func = it.permutations if permutations else it.combinations
        self.pair_inds = np.array(list(it_func(range(self.frag_num), 2)))
        self.frag_inds = np.array([m for m, n in self.pair_inds])

        c3d = geom.coords3d
        frag_c3ds = [c3d[frag] for frag in self.frags]
        self.frag_radii = frag_radii
        if self.frag_radii is None:
            self.frag_radii = [
                get_molecular_radius(frag_c3d, min_offset=radii_offset)
                for frag_c3d in frag_c3ds
            ]
        self.radii_sums = np.array(
            [self.frag_radii[i] + self.frag_radii[j] for i, j in self.pair_inds]
        )

    def get_forces(self, atoms, coords, kappa=None):
        if kappa is None:
            kappa = self.kappa
        c3d = coords.reshape(-1, 3)

        # Break early when only 1 fragment is present
        if len(self.pair_inds) == 0:
            return {"energy": 1, "forces": np.zeros_like(coords)}

        centroids = np.array([c3d[frag].mean(axis=0) for frag in self.frags])
        mm, nn = np.array(self.pair_inds).T
        gdiffs = centroids[mm] - centroids[nn]
        # Add small number to avoid division by zero in 'frag_gradient' calculation
        gnorms = np.linalg.norm(gdiffs, axis=1) + 1e-16
        H = (gnorms < self.radii_sums).astype(int)
        N = H.copy()
        N *= 3 * self.frag_sizes[self.frag_inds]
        N_invs = np.divide(1, N, out=np.zeros_like(N).astype(float), where=N != 0)
        phi = kappa * N_invs * (gnorms - self.radii_sums)
        frag_gradient = (phi * H / gnorms)[:, None] * gdiffs
        gradient = np.zeros_like(c3d)
        # Distribute gradient onto fragments
        for frag_ind, ff in zip(self.frag_inds, frag_gradient):
            gradient[self.frags[frag_ind]] += ff
        forces = -gradient.flatten()

        f3d = forces.reshape(-1, 3)
        f3d -= f3d.mean(axis=0)[None, :]

        return {"energy": 1, "forces": forces}


class PWHardSphere:
    def __init__(
        self,
        geom,
        frags,
        sub_frags,
        kappa=1.0,
    ):
        """Inter-Molecular pairwise Hard-Sphere forces between atoms.

        Hardsphere forces are only applied between certain atoms of given fragments,
        but the whole fragment is moved. Can be used to remove atom inter-molecular
        atom clashes.
        """
        self.frags = frags
        self.sub_frags = sub_frags
        self.kappa = kappa

        self.frag_num = len(self.frags)
        self.frag_sizes = np.array([len(frag) for frag in self.frags])
        self.pair_inds = np.array(list(it.combinations(range(self.frag_num), 2)))
        cov_rads = geom.covalent_radii
        self.pair_atom_inds = list()
        self.pair_cov_radii = list()
        for m, n in self.pair_inds:
            frag_m = self.sub_frags[m]
            frag_n = self.sub_frags[n]
            painds = list()
            pcr = list()
            for fm in frag_m:
                crm = cov_rads[m]
                for fn in frag_n:
                    crn = cov_rads[m]
                    crsum = crm + crn
                    painds.append([fm, fn])
                    pcr.append(crsum)
            self.pair_atom_inds.append(painds)
            self.pair_cov_radii.append(pcr)

    def get_forces(self, atoms, coords, kappa=None):
        if kappa is None:
            kappa = self.kappa
        c3d = coords.reshape(-1, 3)

        # Break early when only 1 fragment is present
        if len(self.pair_inds) == 0:
            return {"energy": 1, "forces": np.zeros_like(coords)}

        forces = np.zeros_like(c3d)
        centroids = np.array([c3d[frag].mean(axis=0) for frag in self.frags])
        N = 1.0
        N_inv = 1 / N
        for (m, n), pai, pcr in zip(
            self.pair_inds, self.pair_atom_inds, self.pair_cov_radii
        ):
            frag_m = self.frags[m]
            frag_n = self.frags[n]
            centr_m = centroids[m]
            centr_n = centroids[n]
            dcentr = centr_m - centr_n
            force_dir = dcentr / np.linalg.norm(dcentr)
            for (am, an), crmn in zip(pai, pcr):
                amn = c3d[am] - c3d[an]
                distmn = np.linalg.norm(amn)
                diff = distmn - crmn
                fact = int(distmn < crmn)
                if not fact:
                    continue
                # Magnitude of applied force
                magn = kappa * N_inv * diff
                force = magn * force_dir
                # Distribute half of the force onto each fragment
                force_2 = force / 2
                # The signs below depend on the difference centr_m - centr_n above
                forces[frag_m] -= force_2
                forces[frag_n] += force_2

        forces -= forces.mean(axis=0)[None, :]
        forces = forces.flatten()

        return {"energy": 1, "forces": forces}
