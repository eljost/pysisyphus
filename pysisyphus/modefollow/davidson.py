#!/usr/bin/env python3

# [1] https://aip.scitation.org/doi/pdf/10.1063/1.1523908
#     Neugebauer, Reiher 2002
# [2] https://reiher.ethz.ch/software/akira.html


from collections import namedtuple

import numpy as np

from pysisyphus.helpers import eigval_to_wavenumber
from pysisyphus.modefollow.NormalMode import NormalMode


def fin_diff(geom, b, step_size):
    m_sqrt = np.sqrt(geom.masses_rep)
    plus = geom.get_energy_and_forces_at(geom.coords + b)["forces"]
    minus = geom.get_energy_and_forces_at(geom.coords - b)["forces"]
    fd = (minus - plus) / (2*step_size) / m_sqrt
    return fd


DavidsonResult = namedtuple(
                    "DavidsonResult",
                    "cur_cycle nus mode_ind",
)


def davidson(geom, q, trial_step_size=0.01, hessian_precon=None, max_cycles=25,
             res_rms_thresh=1e-4):
    print("Using  preconditioner from supplied hessian")
    B_list = list()
    S_list = list()
    msqrt = np.sqrt(geom.masses_rep)

    # Project out translation/rotation
    P = np.eye(geom.cart_coords.size)
    for vec in geom.get_trans_rot_vectors():
        P -= np.outer(vec, vec)
    l_proj = P.dot(q.l_mw) / msqrt
    q = NormalMode(l_proj, geom.masses_rep)

    b_prev = q.l_mw
    for i in range(max_cycles):
        print(f"Cycle {i:02d}")
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
        if hessian_precon is not None:
            # Construct X
            X = np.linalg.inv(hessian_precon - v[mode_ind] * np.eye(hessian_precon.shape[0]))
            b = X.dot(residues[mode_ind])
        else:
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

        if res_rms[mode_ind] < res_rms_thresh:
            print("Converged!")
            break

    result = DavidsonResult(
                cur_cycle=i,
                nus=nus,
                mode_ind=mode_ind,
    )

    return result
