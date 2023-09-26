# [1]   https://doi.org/10.1016/j.cplett.2012.07.003
#       Improved deﬁnition of bond orders for correlated wave functions
#       Mayer, 2012
# [2]   https://doi.org/10.1002/qua.26612
#       Investigating István Mayer's “improved” definitions of bond orders
#       and free valence for correlated singlet-state wave functions
#       Cooper, Ponec, Karadokov, 2021
# [3]   https://doi.org/10.1002/jcc.20494
#       Bond Order and Valence Indices: A Personal Account
#       Mayer, 2006
# [4]   https://doi.org/10.1039/b102094n
#       The Mayer bond order as a tool in inorganic chemistry
#       Bridgeman, Cavigliasso, Ireland, Rothery, 2001


import functools

import numpy as np

from pysisyphus.linalg import matrix_power
from pysisyphus.wavefunction.wavefunction import Wavefunction


@functools.singledispatch
def improved_mayer_bond_order(
    Pa: np.ndarray, Pb: np.ndarray, S: np.ndarray, ao_centers: list, natoms: int
):
    """Improved Mayer bond-order.

    Implements algorithm described in [1].

    Parameters
    ----------
    Pa
        Density matrix of alpha electrons of shape (nao, nao).
    Pb
        Density matrix of beta electrons of shape (nao, nao).
    S
        AO-overlap matrix of shape (nao, nao).
    ao_centers
        List of of length nao, denoting the atom center of each AO.
    natoms
        Number of atoms.

    Returns
    -------
    BOs
        Bond-order matrix.
    """
    D = Pa + Pb
    DS = D @ S

    S_sqrt = matrix_power(S, 0.5)
    S_isqrt = matrix_power(S, -0.5)
    # Eq. (6) in [1]
    uS = 2 * DS - (DS @ DS)
    # Eq. (8) in [1]
    u = S_sqrt @ uS @ S_isqrt
    # Eq. (12) in [1]
    R = S_isqrt @ matrix_power(u, 0.5, thresh=1e-8, strict=False) @ S_isqrt
    RS = R @ S

    BOs = np.zeros((natoms, natoms))
    # Loop over unique basis function pairs
    for i, A in enumerate(ao_centers):
        for j, B in enumerate(ao_centers[i + 1 :], i + 1):
            bo = DS[i, j] * DS[j, i] + RS[i, j] * RS[j, i]
            BOs[A, B] += bo
            BOs[B, A] += bo
    return BOs


@improved_mayer_bond_order.register
def _(wf: Wavefunction):
    """Improved Mayer bond-order.

    Parameters
    ----------
    wf
        Wavefunction object.

    Returns
    -------
    BOs
        Bond-order matrix.
    """

    Pa, Pb = wf.P
    S = wf.S
    ao_centers = list(wf.shells.sph_ao_centers)
    natoms = len(wf.atoms)
    return improved_mayer_bond_order(Pa, Pb, S, ao_centers, natoms)


def report_bond_orders(BOs, thresh=0.1):
    a, b = np.where(BOs > thresh)
    for a_, b_ in zip(a, b):
        if b_ > a_:
            bo = BOs[a_, b_]
            print(f"B({a_: >4d}, {b_: >4d}): {bo:.4f}")
