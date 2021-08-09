# [1] https://doi.org/10.1002/jcc.26495
#     Habershon, 2021

import numpy as np


def get_trans_torque_forces(
    mfrag,
    a_coords3d,
    b_coords3d,
    a_mats,
    b_mats,
    m,
    frags,
    N_inv,
    weight_func=None,
    skip=True,
    kappa=1,
    do_trans=True,
):
    mcoords3d = a_coords3d[mfrag]
    gm = mcoords3d.mean(axis=0)

    if weight_func is None:

        def weight_func(m, n, a, b):
            return 1

    trans_vec = np.zeros(3)
    rot_vec = np.zeros(3)
    for n, nfrag in enumerate(frags):
        if skip and (m == n):
            continue
        amn = a_mats[(m, n)]
        bnm = b_mats[(n, m)]
        for a in amn:
            for b in bnm:
                rd = b_coords3d[b] - a_coords3d[a]
                gd = a_coords3d[a] - gm
                weight = weight_func(m, n, a, b)

                rot_vec += weight * np.cross(rd, gd)
                if do_trans:
                    trans_vec += weight * abs(rd.dot(gd)) * rd / np.linalg.norm(rd)
    trans_vec *= N_inv
    rot_vec *= N_inv
    forces = kappa * (np.cross(-rot_vec, mcoords3d - gm) + trans_vec[None, :])
    return forces


class TransTorque:
    def __init__(
        self,
        frags,
        iter_frags,
        a_mats,
        b_mats,
        weight_func=None,
        skip=True,
        kappa=1.0,
        b_coords3d=None,
        do_trans=True,
    ):
        """Translational and torque forces.
        See A.4. [1], Eqs. (A3) - (A5).
        """
        self.frags = frags
        self.iter_frags = iter_frags
        self.a_mats = a_mats
        self.b_mats = b_mats
        self.weight_func = weight_func
        self.kappa = kappa
        self.skip = skip
        self.b_coords3d = b_coords3d
        self.do_trans = do_trans

        self.set_N_invs()

    def set_N_invs(self):
        Ns = np.zeros(len(self.frags))
        for m, mfrag in enumerate(self.frags):
            for n, _ in enumerate(self.iter_frags):
                if self.skip and (m == n):
                    continue
                amn = self.a_mats[(m, n)]
                bnm = self.b_mats[(n, m)]
                Ns[m] += len(amn) * len(bnm)
            Ns[m] *= 3 * len(mfrag)
        self.N_invs = np.divide(1, Ns, out=np.zeros_like(Ns), where=Ns != 0)

    def get_forces(self, atoms, coords, kappa=None):
        if kappa is None:
            kappa = self.kappa
        c3d = coords.reshape(-1, 3)
        forces = np.zeros_like(c3d)

        b_coords3d = c3d if self.b_coords3d is None else self.b_coords3d

        for m, mfrag in enumerate(self.frags):
            N_inv = self.N_invs[m]
            tt_forces = get_trans_torque_forces(
                mfrag,
                c3d,
                b_coords3d,
                self.a_mats,
                self.b_mats,
                m,
                self.iter_frags,
                N_inv,
                weight_func=self.weight_func,
                skip=self.skip,
                do_trans=self.do_trans,
                kappa=kappa,
            )
            forces[mfrag] = tt_forces

        return {"energy": 1, "forces": forces.flatten()}

    def get_forces_naive(self, atoms, coords, kappa=None):
        if kappa is None:
            kappa = self.kappa

        AR = self.a_mats
        Ns = np.zeros(len(self.frags))
        for m, mfrag in enumerate(self.frags):
            for n, nfrag in enumerate(self.frags):
                if m == n:
                    continue
                Ns[m] += len(AR[(n, m)]) * len(AR[(m, n)])
            Ns[m] *= 3 * len(mfrag)

        c3d = coords.reshape(-1, 3).copy()
        f3d = np.zeros_like(c3d)
        vts = np.zeros((len(self.frags), 3))
        vrs = np.zeros((len(self.frags), 3))
        for m, mfrag in enumerate(self.frags):
            gm = c3d[mfrag].mean(axis=0)
            c3dm = c3d[mfrag]
            for n, nfrag in enumerate(self.frags):
                if m == n:
                    continue
                for a in AR[(m, n)]:
                    for b in AR[(n, m)]:
                        rdiff = c3d[b] - c3d[a]
                        rdiffn = np.linalg.norm(rdiff)
                        rg = c3d[a] - gm
                        quot = rdiff / rdiffn
                        dot = rdiff.dot(rg)
                        vts[m] += abs(dot) * quot
                        vrs[m] += np.cross(rdiff, rg)
            vts[m] /= Ns[m]
            vrs[m] /= Ns[m]
            if not self.do_trans:
                vts[m] = 0.0
            f3d[mfrag] = kappa * (-np.cross(vrs[m], c3dm-gm[None, :]) + vts[m])

        forces = f3d.flatten()
        return {"energy": 1, "forces": forces}

    # def get_forces(self, atoms, coords, kappa=None):
        # return self.get_forces_naive(atoms, coords, kappa=kappa)
