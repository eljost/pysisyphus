# [1] https://doi.org/10.1002/jcc.540140615
#     Comparison of the Boys and Pipek–Mezey localizations in the local
#     correlation approach and automatic virtual basis selection
#     Boughton, Pulay, 1993
# [2] https://doi.org/10.1063/1.2360264
#     Fast noniterative orbital localization for large molecules
#     Aquilante, Pedersen, 2006


from functools import singledispatch
import itertools as it
from typing import List

import numpy as np
from numpy.typing import NDArray

from pysisyphus.wavefunction import logger, Wavefunction
from pysisyphus.linalg import pivoted_cholesky


@singledispatch
def cholesky(C: NDArray[float]):
    """Localization via pivoted Cholesky factorization.

    See [2].

    Parameters
    ----------
    C
        Matrix of molecular orbital coefficients to be localized.
        Shape is (naos, nmos).

    Returns
    -------
    C_loc
        Localized molecular orbital coefficient matrix of shape (naos, nmos).
    """
    naos, nmos = C.shape

    # We can't use scipy's builtin Cholesky-factorization, as it does not
    # support positive semi-definite matrices.
    L, piv, _ = pivoted_cholesky(C @ C.T)
    # Restore original ordering via permutation matrix.
    P = np.zeros((naos, naos))
    P[piv, np.arange(naos)] = 1
    C_loc = P @ L[:, :nmos]
    return C_loc


@cholesky.register
def _(wf: Wavefunction):
    """Currently localizes only C_(α,occ) MOs."""
    Cao, _ = wf.C_occ
    return cholesky(Cao)


def rot_inplace(mat, rad, i, j):
    """Inplace rotation of matrix columns A[:, i] and A[:, j] by 'rad' radians."""
    cos = np.cos(rad)
    sin = np.sin(rad)
    i_rot = cos * mat[:, i] + sin * mat[:, j]
    j_rot = -sin * mat[:, i] + cos * mat[:, j]
    mat[:, i] = i_rot
    mat[:, j] = j_rot


@singledispatch
def pipek_mezey(
    C: NDArray[float],
    S: NDArray[float],
    ao_center_map: dict[int, List[int]],
    max_cycles: int = 100,
    dP_thresh: float = 1e-8,
):
    """Pipek-Mezey localization using Mulliken population analysis.

    Python adaption of code found in orbloc.f90 of Multiwfn 3.8.
    For now, only Mulliken population and localization exponent 2 is supported.

    Parameters
    ----------
    C
        Matrix of molecular orbital coefficients to be localized.
        Shape is (naos, nmos).
    S
        Overlap matrix of shape (naos x naos).
    ao_center_map
        Mapping between atom indices and AOs, centered at the respective atom.
    max_cycles
        Maximum number of macro cycles to achieve successful localization.
    dp_thresh
        Indicate convergence when change in cost function is equal or below
        this threshold.

    Returns
    -------
    C_loc
        Localized molecular orbital coefficient matrix of shape (naos, nmos).
    """

    C = C.copy()
    _, nmos = C.shape
    centers = list(ao_center_map.keys())

    pi_quart = np.pi / 4
    SC = S @ C
    P_prev = 0.0
    logger.info(f"Pipek-Mezey localization")
    for i in range(max_cycles):
        # Loop over pairs of MO indices and do 2x2 rotations.
        for j, k in it.combinations(range(nmos), 2):
            Q = SC * C
            A = 0.0
            B = 0.0
            # Eq. (9) in [1]
            for center in centers:
                ao_inds = ao_center_map[center]
                Qjk = (
                    C[ao_inds, j] * SC[ao_inds, k] + C[ao_inds, k] * SC[ao_inds, j]
                ).sum() / 2
                Qjj = Q[ao_inds, j].sum()
                Qkk = Q[ao_inds, k].sum()
                A += (Qjk ** 2) - ((Qjj - Qkk) ** 2) / 4
                B += Qjk * (Qjj - Qkk)

            if (A ** 2 + B ** 2) <= 1e-12:
                continue

            gamma = np.sign(B) * np.arccos(-A / np.sqrt(A ** 2 + B ** 2)) / 4
            assert -pi_quart <= gamma <= pi_quart
            # 2x2 MO rotations
            rot_inplace(C, gamma, j, k)
            # The MO rotations invalidate the SC matrix product. We update it too.
            rot_inplace(SC, gamma, j, k)
        # Outside of loop over orbital pairs

        # We wan't to maximize P and we monitor the progress. Eq. (8) in [1].
        P = 0.0
        for l in range(nmos):
            for center in centers:
                ao_inds = ao_center_map[center]
                Q = (C[ao_inds, l][:, None] * C[:, l] * S[ao_inds, :]).sum()
                P += Q ** 2
        dP = P - P_prev
        logger.info(f"{i:03d}: {P=: >12.8f} {dP=: >12.8f}")
        if dP <= dP_thresh:
            logger.info("Converged after {i+1} cycles.")
            break
        P_prev = P
    # Outside macro cycles
    return C


@pipek_mezey.register
def _(wf: Wavefunction):
    """Currently localizes only C_(α,occ) MOs."""
    S = wf.S
    Cao, _ = wf.C_occ
    C_chol_loc = cholesky(Cao)
    return pipek_mezey(C_chol_loc, S, wf.ao_center_map)
