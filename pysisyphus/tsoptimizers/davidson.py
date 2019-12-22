#!/usr/bin/env python3

# [1] https://aip.scitation.org/doi/pdf/10.1063/1.1523908
#     Neugebauer, Reiher 2002
# [2] https://reiher.ethz.ch/software/akira.html


from math import sqrt
import sys

import numpy as np

from pysisyphus.calculators.Turbomole import Turbomole
from pysisyphus.helpers import geom_from_xyz_file, eigval_to_wavenumber


class NormalMode:
    """See http://gaussian.com/vib/"""

    def __init__(self, l, masses):
        """NormalMode class.

        Cartesian displacements are normalized to 1.

        Parameters
        ----------
        l : np.array
            Cartesian, non-mass-weighted displacements.
        masses : np.array
            Atomic masses.
        """

        self.l = np.array(l)
        self.l /= np.linalg.norm(l)
        self.masses = masses
        assert self.l.shape == self.masses.shape

    @property
    def red_mass(self):
        return 1 / np.sum(np.square(self.l_mw) / self.masses)

    def mw_norm_for_norm(self, norm=0.01):
        return norm * sqrt(self.red_mass)

    @property
    def l_mw(self):
        l_mw = self.l * np.sqrt(self.masses)
        return l_mw / np.linalg.norm(l_mw)


def fin_diff(geom, b, step_size):
    m_sqrt = np.sqrt(geom.masses_rep)
    plus = geom.get_energy_and_forces_at(geom.coords + b)["forces"]
    minus = geom.get_energy_and_forces_at(geom.coords - b)["forces"]
    fd = (minus - plus) / (2*step_size) / m_sqrt
    return fd


def davidson(geom, q, trial_step_size=0.01):
    B_list = list()
    S_list = list()
    msqrt = np.sqrt(geom.masses_rep)

    P = np.eye(geom.cart_coords.size)
    for vec in geom.get_trans_rot_vectors():
        P -= np.outer(vec, vec)

    # Project out translation/rotation
    l_proj = P.dot(q.l_mw) / msqrt
    q = NormalMode(l_proj, geom.masses_rep)
    b_prev = q.l_mw

    for i in range(10):
        b = q.l_mw
        B_list.append(b)
        B = np.stack(B_list, axis=1)

        # Overlaps of basis vectors in B
        # B_ovlp = np.einsum("ij,kj->ik", B_list, B_list)
        # print("Basis vector overlaps")
        # print(B_ovlp)

        # Estimate s (sigma) from finite differences 

        # Get step size in mass-weighted coordinates that results 
        # in the desired 'trial_step_size' in not-mass-weighted coordinates.
        mw_step_size = q.mw_norm_for_norm(trial_step_size)
        # Actual step in non-mass-weighted coordinates
        step = trial_step_size * q.l
        s = fin_diff(geom, step, mw_step_size)
        S_list.append(s)
        S = np.stack(S_list, axis=1)

        # Calculate and symmetrize approximate hessian
        Hm_ = B.T.dot(S)
        Hm_ = (Hm_ + Hm_.T) / 2

        # Diagonalization
        v, w = np.linalg.eigh(Hm_)

        # i-th approximation to exact eigenvector
        approx_modes = (w*B[:,:,None]).sum(axis=1).T

        # Calculate overlaps between previous root and the new approximate
        # normal modes for root following.
        mode_overlaps = (approx_modes * b_prev).sum(axis=1)
        mode_ind = np.abs(mode_overlaps).argmax()
        print(f"\tFollowing mode {mode_ind}")

        residues = list()
        for s in range(i+1):
            residues.append(
                (w[:,s] * (S-v[s]*B)).sum(axis=1)
            )
        residues = np.array(residues)

        b_prev = approx_modes[mode_ind]

        # Construct new basis vector from residuum of selected mode
        b = residues[mode_ind]
        # Project out translation and rotation from new mode guess
        b = P.dot(b)
        # Orthogonalize new mode against current basis vectors
        rows, cols = B.shape
        B_ = np.zeros((rows, cols+1))
        B_[:,:cols] = B
        B_[:,-1] = b
        b, _ = np.linalg.qr(B_)

        # New NormalMode from non-mass-weighted displacements
        q = NormalMode(b[:,-1] /msqrt, geom.masses_rep)

        # Calculate wavenumbers
        nus = eigval_to_wavenumber(v)

        # Check convergence criteria
        max_res = np.abs(residues).max(axis=1)
        res_rms = np.sqrt(np.mean(residues**2, axis=1))

        # Print progress
        print("\t #  |      wavelength       |  rms       |   max")
        for j, (nu, rms, mr) in enumerate(zip(nus, res_rms, max_res)):
            sel_str = "*" if (i == mode_ind) else " "
            print(f"\t{j:02d}{sel_str} | {nu:> 16.2f} cm⁻¹ | {rms:.8f} | {mr:.8f}")
        print()

        sys.stdout.flush()
        if res_rms[mode_ind] < 1e-4:
            print("Converged!")
            break

    assert (nus[mode_ind] - 1748.21916159) < 1e-4


def test():
    geom = geom_from_xyz_file("/scratch/programme/pysisyphus/pysisyphus/tsoptimizers/coords.xyz")
    control_path = "/scratch/programme/pysisyphus/pysisyphus/tsoptimizers/acet_tm"
    geom.set_calculator(Turbomole(control_path=control_path))
    l = np.zeros_like(geom.coords).reshape(-1, 3)
    l[0][2] = 0.8
    l[1][2] = -0.6
    q = NormalMode(l.flatten(), geom.masses_rep)
    
    davidson(geom, q)


if __name__ == "__main__":
    test()
