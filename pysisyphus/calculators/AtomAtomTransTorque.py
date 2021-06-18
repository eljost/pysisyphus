# [1] https://doi.org/10.1002/jcc.26495
#     Habershon, 2021

import itertools as it

import numpy as np

from pysisyphus.elem_data import COVALENT_RADII as CR


class AtomAtomTransTorque:
    def __init__(
        self,
        geom,
        frags,
        A_mats,
        kappa=2.0,
    ):
        """Atom-atom translational and torque forces.

        See A.5. [1], Eq. (A6).
        """
        self.geom = geom
        self.frags = frags
        self.A_mats = A_mats
        self.kappa = kappa

        self.frag_num = len(self.frags)
        self.pair_inds = list(it.permutations(range(self.frag_num), 2))
        self.frag_sizes_sq = [len(self.frags[m]) ** 2 for m, _ in self.pair_inds]

        frag_atoms = [[self.geom.atoms[a].lower() for a in frag] for frag in self.frags]
        frag_cov_rads = [np.array([CR[fa.lower()] for fa in fas]) for fas in frag_atoms]
        self.avg_cov_radii = list()
        for m, n in self.pair_inds:
            m_cov_rads = frag_cov_rads[m]
            n_cov_rads = frag_cov_rads[n]
            mn_avg_cov_radii = (m_cov_rads[:, None] + n_cov_rads[None, :]) / 2
            self.avg_cov_radii.append(mn_avg_cov_radii)

        def phi_func(A, a):
            """See (A6) in [1]."""
            return 1.5 if a in A else 2.0

        self.phis = list()
        for m, n in self.pair_inds:
            key = (m, n)
            A = self.A_mats[key]
            self.phis.append(np.array([phi_func(A, a) for a in self.frags[m]]))

    # def get_forces_naive(self, atoms, coords):
        # def p(s):
            # print(f"REF: {s}")

        # c3d = coords.reshape(-1, 3)
        # gs = [c3d[frag].mean(axis=0) for frag in self.frags]

        # zt = np.zeros((self.frag_num, 3))
        # zr = np.zeros((self.frag_num, 3))
        # for (m, n), phis in zip(self.pair_inds, self.phis):
            # p(f"m={m}, n={n}")
            # mfrag = self.frags[m]
            # nfrag = self.frags[n]
            # A = self.A_mats[(m, n)]
            # g = gs[m]
            # N = 0
            # for a in mfrag:
                # ra = c3d[a]
                # ramg = ra - g
                # for b in nfrag:
                    # rb = c3d[b]
                    # rmr = rb - ra
                    # a_atm = atoms[a].lower()
                    # b_atm = atoms[b].lower()
                    # cab = (CR[a_atm] + CR[b_atm]) / 2
                    # phi = 1.5 if a in A else 2
                    # nrmr = np.linalg.norm(rmr)
                    # H = 1 if nrmr < (phi * cab) else 0
                    # y = cab * rmr / np.linalg.norm(rmr) - rmr
                    # p(f"y_(a={a},b={b})={y}")
                    # ny = np.linalg.norm(y)
                    # zt[m] += H * y.dot(ramg) * y / ny
                    # p(f"dot={y.dot(ramg)}")
                    # zr[m] += H * np.cross(y, ramg)
                    # p(f"zt={zt}")
                    # p(f"zr={zr}")
                    # N += H
            # N *= self.frag_sizes_sq[m]
            # p(f"N={N}")
            # if N == 0:
                # N_inv = 0
            # else:
                # N_inv = 1 / N
            # zt[m] *= N_inv
            # zr[m] *= N_inv

        # forces = np.zeros_like(c3d)
        # for m, mfrag in enumerate(self.frags):
            # forces[mfrag] = np.cross(-zr[m], c3d[mfrag] - gs[m]) + zt[m]
        # forces *= self.kappa
        # # return zt, zr
        # print("REF")
        # print(forces)
        # return {"energy": 1, "forces": forces.flatten()}

    def get_forces(self, atoms, coords):
        c3d = coords.reshape(-1, 3)

        if len(self.pair_inds) == 0:
            return {"energy": 1, "forces": np.zeros_like(coords)}

        frag_coords = [c3d[frag] for frag in self.frags]
        frag_centroids = np.array([c3d[frag].mean(axis=0) for frag in self.frags])
        frag_centered = [
            frag_coords[m] - frag_centroids[m] for m in range(self.frag_num)
        ]
        Hs = np.zeros(self.frag_num)

        zt = np.zeros((self.frag_num, 3))
        zr = np.zeros((self.frag_num, 3))
        for (m, n), phi, avg_cov_radii in zip(
            self.pair_inds, self.phis, self.avg_cov_radii
        ):
            mfrag = self.frags[m]
            mcoords3d = c3d[mfrag]
            ncoords3d = c3d[self.frags[n]]
            coord_diffs = ncoords3d[None, :] - mcoords3d[:, None]
            norms = np.linalg.norm(coord_diffs, axis=2)
            H = norms < phi[:, None] * avg_cov_radii
            H_sum = H.sum()
            if H_sum == 0:
                continue
            Hs[m] += H_sum
            y = (avg_cov_radii / norms - 1)[:, :, None] * coord_diffs
            ynorms = np.linalg.norm(y, axis=2)
            mcoords3d_centered = frag_centered[m]
            dot = np.abs(np.einsum("ijk,ik->ij", y, mcoords3d_centered))
            zt[m] += (H[:, :, None] * dot[:, :, None] * y / ynorms[:, :, None]).sum(
                axis=(0, 1)
            )
            zr[m] += (H[:, :, None] * np.cross(y, mcoords3d_centered[:, None, :])).sum(
                axis=(0, 1)
            )
        N = self.frag_sizes_sq * Hs
        N_invs = np.divide(1, N, out=np.zeros_like(N), where=N != 0)
        zt *= N_invs[:, None]
        zr *= N_invs[:, None]

        forces = np.zeros_like(c3d)
        for m, mfrag in enumerate(self.frags):
            forces[mfrag] = np.cross(-zr[m], frag_centered[m]) + zt[m]
        forces *= self.kappa

        return {"energy": 1, "forces": forces.flatten()}
