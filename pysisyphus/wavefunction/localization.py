# [1] https://doi.org/10.1002/jcc.540140615
#     Comparison of the Boys and Pipek–Mezey localizations in the local
#     correlation approach and automatic virtual basis selection
#     Boughton, Pulay, 1993
# [2] https://doi.org/10.1063/1.2360264
#     Fast noniterative orbital localization for large molecules
#     Aquilante, Pedersen, 2006
# [3] https://doi.org/10.1063/1.3042233
#     Constructing diabatic states from adiabatic states: Extending
#     generalized Mulliken–Hush to multiple charge centers with Boys localization
#     Subotnik, Yeganeh, Cave, Ratner, 2008
# [4] https://doi.org/10.1063/1.1681683
#     Localized molecular orbitals for polyatomic molecules.
#     I. A comparison of the Edmiston‐Ruedenberg and Boys localization methods
#     Kleier, Halgren, Hall Jr., Lipscomb, 1974
# [5] https://doi.org/10.1063/1.4894472
#     Diabatization based on the dipole and quadrupole: The DQ method
#     Hoyer, Xu, Ma, Gagliardi, Truhlar, 2014
# [6] https://pubs.acs.org/doi/pdf/10.1021/acs.jctc.2c00261
#     Implementation of Occupied and Virtual Edmiston−Ruedenberg
#     Orbitals Using Cholesky Decomposed Integrals
# [7] https://doi.org/10.1063/1.1790971
#     An efficient method for calculating maxima of homogeneous
#     functions of orthogonal matrices: Applications to localized
#     occupied orbitals


from dataclasses import dataclass
from functools import singledispatch
import itertools as it
from typing import Callable, Optional
import warnings

import numpy as np
from numpy.typing import NDArray

from pysisyphus.helpers_pure import rms
from pysisyphus.linalg import matrix_power, pivoted_cholesky
from pysisyphus.wavefunction import logger, Wavefunction
from pysisyphus.wavefunction.DIIS import DIIS


PI_QUART = np.pi / 4


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
    """Inplace rotation of matrix columns mat[:, i] and mat[:, j] by 'rad' radians."""
    cos = np.cos(rad)
    sin = np.sin(rad)
    i_rot = cos * mat[:, i] + sin * mat[:, j]
    j_rot = -sin * mat[:, i] + cos * mat[:, j]
    mat[:, i] = i_rot
    mat[:, j] = j_rot


@dataclass
class JacobiSweepResult:
    is_converged: bool
    cur_cycle: int
    C: NDArray[float]
    P: float


def jacobi_sweeps(
    C: NDArray[float],
    cost_func: Callable,
    ab_func: Callable,
    callback: Optional[Callable] = None,
    max_cycles: int = 100,
    dP_thresh: float = 1e-8,
) -> JacobiSweepResult:
    """Wrapper for 2x2 Jacobi-sweeps as used in localization/diabatization.

    Parameters
    ----------
    C
        MO coefficient matrix (shape naos x nmos) or rotation matrix (nstates x nstates).
    cost_func
        Function to be maximized/minimized.
    ab_func
        Function that returns A & B values, used to calculate the angle
        for the 2x2 rotation.
    callback
        Function that is called after the 2x2 rotation took place. It takes three arguments:
        (gamma, j, k).
    max_cycles
        Maximum number of macro cycles.
    dP_thresh
        Indicate convergence when change in cost function is equal or below
        this threshold.

    Returns
    -------
    C_loc
        Localized molecular orbital coefficient matrix of shape (naos, nmos).
    """

    assert max_cycles > 0
    C = C.copy()
    _, nmos = C.shape

    if callback is None:
        callback = lambda *args: None

    P_prev = 0.0

    logger.info(f"Starting Jacobi sweeps.")
    for i in range(max_cycles):
        # Loop over pairs of MO indices and do 2x2 rotations.
        for j, k in it.combinations(range(nmos), 2):
            A, B = ab_func(j, k, C)

            if (A**2 + B**2) <= 1e-12:
                continue

            gamma = np.sign(B) * np.arccos(-A / np.sqrt(A**2 + B**2)) / 4
            assert -PI_QUART <= gamma <= PI_QUART
            # 2x2 MO rotations
            rot_inplace(C, gamma, j, k)
            callback(gamma, j, k)
        # Outside of loop over orbital pairs

        # Calculate the target cost function we want to maximize/minimize.
        P = cost_func(C)
        dP = P - P_prev
        logger.info(f"{i:03d}: {P=: >12.8f} {dP=: >12.8f}")
        if is_converged := (dP <= dP_thresh):
            logger.info(f"Jacobi sweeps converged in {i+1} cycles.")
            break
        P_prev = P
    # Outside macro cycles

    result = JacobiSweepResult(
        is_converged=is_converged,
        cur_cycle=i,
        C=C,
        P=P,
    )
    return result


@singledispatch
def pipek_mezey(C, S, ao_center_map, **kwargs) -> JacobiSweepResult:
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

    Returns
    -------
    C_loc
        Localized molecular orbital coefficient matrix of shape (naos, nmos).
    """
    _, nmos = C.shape
    centers = list(ao_center_map.keys())

    SC = S @ C

    def ab_func(j, k, C):
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
            A += (Qjk**2) - ((Qjj - Qkk) ** 2) / 4
            B += Qjk * (Qjj - Qkk)
        return A, B

    def callback(gamma, j, k):
        # The MO rotations invalidate the SC matrix product. We update it too.
        rot_inplace(SC, gamma, j, k)

    def cost_func(C):
        # We wan't to maximize P and we monitor the progress. Eq. (8) in [1].
        P = 0.0
        for l in range(nmos):
            for center in centers:
                ao_inds = ao_center_map[center]
                Q = (C[ao_inds, l][:, None] * C[:, l] * S[ao_inds, :]).sum()
                P += Q**2
        return P

    logger.info("Pipek-Mezey localization")
    return jacobi_sweeps(C, cost_func, ab_func, callback, **kwargs)


@pipek_mezey.register
def _(wf: Wavefunction) -> JacobiSweepResult:
    """Currently localizes only C_(α,occ) MOs."""
    S = wf.S
    Cao, _ = wf.C_occ
    C_chol_loc = cholesky(Cao)
    return pipek_mezey(C_chol_loc, S, wf.ao_center_map)


def get_fb_contract(moments_ints):
    def contract(mo_j, mo_k):
        return np.einsum(
            "xkl,k,l->x",
            moments_ints,
            mo_j,
            mo_k,
            optimize=["einsum_path", (0, 1), (0, 1)],
        )

    return contract


def get_fb_ab_func(moments_ints):
    contract = get_fb_contract(moments_ints)

    def ab_func(j, k, C):
        mo_j = C[:, j]
        mo_k = C[:, k]
        jrj = contract(mo_j, mo_j)
        jrk = contract(mo_j, mo_k)
        krk = contract(mo_k, mo_k)
        A = (jrk**2).sum() - ((jrj - krk) ** 2).sum() / 4  # Eq. (9) in [4]
        B = ((jrj - krk) * jrk).sum()  # Eq. (10) in [4]
        return A, B

    return ab_func


def get_fb_cost_func(moments_ints):
    def cost_func(C):
        vals = np.einsum(
            "xkl,ki,li->ix",
            moments_ints,
            C,
            C,
            optimize=["einsum_path", (0, 1), (0, 1)],
        )
        val = ((vals[:, None, :] - vals) ** 2).sum()
        return val

    return cost_func


@singledispatch
def foster_boys(
    C: NDArray[float], dip_ints: NDArray[float], **kwargs
) -> JacobiSweepResult:
    """Foster-Boys localization.

     nMO   nMO
     ___   ___
     ╲     ╲                         2
      ╲     ╲   |                   |
      ╱     ╱   | <i|R|i> - <j|R|j> |
     ╱     ╱    |                   |
     ‾‾‾   ‾‾‾
    i = 1 j = 1

    or similarily (see eq. (6) in [4] or the appendix of [5])

     nMO
     ___
     ╲                2
      ╲    |         |
      ╱    | <i|R|i> |
     ╱     |         |
     ‾‾‾
    i = 1

    Parameters
    ----------
    dip_ints
        Dipole moment integral matrices with shape (3, naos, naos).
    C
        Matrix of molecular orbital coefficients to be localized.
        Shape is (naos, nmos).
    kwargs
        Additional keyword arguments that are passed jacobi_sweeps.

    Returns
    -------
    C_loc
        Localized molecular orbital coefficient matrix of shape (naos, nmos).
    """

    ab_func = get_fb_ab_func(dip_ints)
    cost_func = get_fb_cost_func(dip_ints)
    logger.info("Foster-Boys localization")
    return jacobi_sweeps(C, cost_func, ab_func, **kwargs)


@foster_boys.register
def _(wf: Wavefunction) -> JacobiSweepResult:
    """Currently localizes only C_(α,occ) MOs."""
    Cao, _ = wf.C_occ
    dip_ints = wf.dipole_ints()
    C_chol_loc = cholesky(Cao)
    return foster_boys(C_chol_loc, dip_ints)


def edmiston_ruedenberg_cost_func(L):
    """Edmiston-Ruedenberg cost function.

    Eq. (2) in [7]."""
    return np.einsum("Jpp,Jpp->", L, L, optimize="greedy")


def edmiston_ruedenberg_grad(L):
    """Edmiston-Ruedenberg gradient.

    Eq. (5) in [7]."""
    g1 = np.einsum("Jpq,Jqq->pq", L, L, optimize="greedy")
    g2 = np.einsum("Jqp,Jpp->pq", L, L, optimize="greedy")
    return -4 * (g1 - g2)


def edmiston_ruedenberg(C, Lao, Sao, diis_cycles=5, max_cycles=100):
    """Edmiston-Ruedenberg localization using DF integrals and DIIS.

    Parameters
    ----------
    C
        MO-coefficients of shape (naos, nmos); this means MOs should
        be in columns.

    Lao
        Density fitting tensor obtained from contracting the 2e3c-matrix
        with 2c2e**-0.5 in the AO basis. Must have shape (naux, nao, nao).
        Same as in
    Sao
        Overlap matrix of shape (nao, nao) in the AO basis.
    diis_cycles
        Integer >= 0. If > 0, DIIS is employed to accelerate convergence.
    max_cycles
        Positive integer; maximum number of localization cycles.

    Returns
    -------
    Crot
        Edmiston-Ruedenberg localized orbitals.
    """
    nmos = C.shape[1]
    # We can't keep more error vectors than MOs
    if nmos < diis_cycles:
        warnings.warn(
            f"Can keep at most {nmos} DIIS error vectors, but "
            f"{diis_cycles} were requested!."
        )
    diis_cycles = min(diis_cycles, nmos)
    if diis_cycles == 1:
        diis_cycles = 0
    started_to_store_diis = False
    specs = {
        "err_vecs": np.zeros((diis_cycles, nmos**2)),
        "R_mats": np.zeros((diis_cycles, nmos, nmos)),
        "D_mats": np.zeros((diis_cycles, nmos, nmos)),
    }
    diis = DIIS(specs)

    D = np.eye(nmos)
    C0 = C.copy()
    Crot = C.copy()
    for i in range(max_cycles):
        # Transform cholesky integrals to MO basis, step (2) in [7].
        Lpq = np.einsum("Luv,up,vq->Lpq", Lao, Crot, Crot, optimize="greedy")
        f = edmiston_ruedenberg_cost_func(Lpq)
        g = edmiston_ruedenberg_grad(Lpq)
        gnorm = np.linalg.norm(g)
        grms = rms(g)
        # Construct transformation, step (3) in [7].
        # Rji = (Xj Xi Xi Xi)
        R = np.einsum("Jji,Jii->ji", Lpq, Lpq, optimize="greedy")
        # DIIS error; lack of symmetry in R.
        err = (R - R.T).flatten()
        errrms = rms(err)
        print(
            f"Cycle {i:02d} f={f:.6f}, |g|={np.linalg.norm(g):.4f}, rms(g)={grms:>8.4e} "
            f"rms(err)={errrms:>8.4e}"
        )
        if converged := grms <= 1e-4:
            print("Edmiston-Ruedenberg localization converged!")
            break

        U = R @ matrix_power(R.T @ R, -0.5)
        D = D @ U
        if diis_cycles and (started_to_store_diis or (gnorm <= 2.0e-1)):
            diis.store(
                {
                    "err_vecs": err,
                    "R_mats": R,
                    "D_mats": D,
                }
            )
            started_to_store_diis = True
        if (diis_coeffs := diis.get_coeffs()) is not None:
            D_diis = np.einsum("i,ijk->jk", diis_coeffs, diis.get("D_mats"))
            C_diis = C0 @ D_diis
            # DIIS-2 algorithm in [7]
            R_diis = np.einsum("i,ijk->jk", diis_coeffs, diis.get("R_mats"))
            S_diis = C_diis.T @ Sao @ C_diis
            S_diis_inv = matrix_power(S_diis, -1.0)
            # Generalized eta step
            V = S_diis_inv @ R_diis @ matrix_power(R_diis.T @ S_diis_inv @ R_diis, -0.5)
            # Recalculate D from DIIS results
            D = D_diis @ V
        # Update Crot, step (4) in [7]. This either uses the D matrix obtained from DIIS
        # or the one previously calculated.
        Crot = C0 @ D
    # TODO: return JacobiSweepResult?
    return Crot, f
