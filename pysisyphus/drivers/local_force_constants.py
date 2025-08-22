# [1] https://doi.org/10.1063/1.4747339
#     Relating normal vibrational modes to local vibrational
#     modes with the help of an adiabatic connection scheme
#     Zou, Kalescky, Kraka, Cremer, 2012
# [2] https://doi.org/10.1016/j.cplett.2020.137337
#     Local vibrational force constants – From the assessment
#     of empirical force constants to the description of bonding
#     in large systems
#     Zou, Tao, Freindorf, Cremer, Kraka, 2020
# [3] https://pubs.acs.org/doi/epdf/10.1021/acs.jctc.8b01279
#     In Situ Measure of Intrinsic Bond Strength
#     in Crystalline Structures: Local Vibrational Mode Theory for
#     Periodic Systems
#     Tao, Zou, Sethio, Verma, Qiu, Tian, Cremer, Kraka, 2019
# [4] https://doi.org/10.1021/acs.jctc.1c01269
#     LModeA-nano: A PyMOL Plugin for Calculating Bond Strength
#     in Solids, Surfaces, and Molecules via Local Vibrational Mode Analysis
#     Tao, Zou, Nanayakkara, Kraka, 2022
#
# I actually never read [3] and [4], but they seem to give a good overview.


import functools
import warnings

import numpy as np


from pysisyphus.constants import AU2MDYNEPERANG
from pysisyphus.Geometry import Geometry


def compliance_mat(hessian, B):
    """Compliance matrix from Cartesian Hessian and Wilson-B matrix.

    See eqs. (6), (10) and (14) in [1]:
    Eq. (6):  K = L^T f^x L
    Eq. (10): D = B L
    Eq. (14): Γ_q = (F_q)⁻¹ = D K⁻¹ D^T

    (10) in (14):
    Γ_q = B L K⁻¹ (B L)^T
        = B L K⁻¹ L^T B^T

    From (6):
    K = L^T f^x L
    K L^T = L^T f^x
    L K L^T = f^x
    L K⁻¹ L^T = f^x⁻¹

    So: Γ_q = B f^x⁻¹

    """
    hessian_inv = np.linalg.pinv(hessian)
    return B @ hessian_inv @ B.T


def get_force_constants_from_complice_mat(hessian, B):
    """
    Parameters
    ----------
    hessian
        Cartesian Hessian matrix.
    B
        Wilson-B-matrix.

    Returns
    -------
    force_constants
        Local force constants in atomic units (mass/time**2).
    gamma_q
        Compliance matrix.
    """
    gamma_q = compliance_mat(hessian, B)
    force_constants = 1 / np.diag(gamma_q)
    return force_constants, gamma_q


def get_local_force_constants(hessian, B, L):
    """Local force constants and local normal modes.

    Parameters
    ----------
    hessian
        Cartesian Hessian matrix, shape (ncart, ncart).
    B
        Wilson-B-matrix, shape (ninternsl, ncart).
    L
        Renormalized Cartesian displacment vectors in columns.
        Obtained from unweighting normal modes and renormalizing.
        Shape (ncart, nnmodes).

    Returns
    -------
    force_constants
        Local force constants in atomic units (mass/time**2).
    local_modes
        Local vibrational modes in the basis of the normal modes.
        Local modes are given in columns, so the shape of local_modes is
        (nmodes, ninternal).
    """

    ninternal = B.shape[0]
    nnmodes = L.shape[1]
    # Eq. (10) in [1]; D has shape (ninternal, nmodes)
    D = B @ L
    # Eq. (6) in [1]
    K = L.T @ hessian @ L
    # Diagonal of inverted K matrix.
    diag_K_inv = 1 / np.diag(K)
    force_constants = np.zeros(ninternal)
    local_modes = np.zeros((nnmodes, ninternal))
    for i, d_row in enumerate(D):
        # As K_inv is a diagonal matrix we can avoid the matrix multiplications.
        K_inv_dT = diag_K_inv * d_row
        dK_inv_dT = d_row.dot(K_inv_dT)

        # Eq. (15) in [2]
        force_constants[i] = 1 / dK_inv_dT
        # Eq. (12) in [2]
        local_modes[:, i] = K_inv_dT * force_constants[i]
    return force_constants, local_modes


def local_mode_overlaps(hessian, L, local_modes):
    """
    Parameters
    ----------
    hessian
        Cartesian Hessian matrix.
    L
        Renormalized Cartesian displacment vectors. Obtained from unweighting
        normal modes and renormalizing.
    local_modes
        Local vibrational modes.

    Returns
    -------
    S
        Overlap matrix of shape (n_local_modes, n_normal_modes).
    C
        Decomposition matrix of shape (n_local_modes, n_normal_modes). Every
        column belongs to one normal mode and contains the contributions of
        each local normal mode to this normal mode. Each column sums to 1.0.
    """
    # Eq. (28) in [2]
    cart_local_modes = L @ local_modes
    # Matrix in numerator in eq (27) in [2]
    anlm = np.einsum("in,ij,jm->nm", cart_local_modes, hessian, L, optimize="greedy")
    # Vectors in denominator of eq (27) in [2]
    anan = np.einsum(
        "in,ij,jn->n", cart_local_modes, hessian, cart_local_modes, optimize="greedy"
    )
    lmlm = np.einsum("im,ij,jm->m", L, hessian, L, optimize="greedy")

    # Eq. (27) in [2]
    # S has shape (nlocalmodes, nnormalmodes)
    S = anlm**2 / anan[:, None] / lmlm[None, :]

    # Eq. (30) in [2]
    C = S / S.sum(axis=0)[None, :]
    return S, C


@functools.singledispatch
def stretch_constants(hessian: np.ndarray, L, typed_prims, B):
    force_constants, local_modes = get_local_force_constants(hessian, B, L)
    results = dict()
    for i, (tp, fc) in enumerate(zip(typed_prims, force_constants)):
        fc_cgs = fc * AU2MDYNEPERANG
        # print(f"{i:03d}: {tp}, k={fc_cgs:8.3f} mdyn/Å")
        results[tp] = fc
    return results


@stretch_constants.register
def _(geom: Geometry, proj_gradient: bool = True, **kwargs):
    hessian = geom.cart_hessian
    # Please note, that we also project out the remaining gradient. This will only be
    # done at non-stationary points with non-vanishing gradient.
    nus, *_, L = geom.get_normal_modes(proj_gradient=proj_gradient)
    neg_mask = nus <= 0.0
    if neg_mask.sum() > 0:
        neg_nus = nus[neg_mask]
        warnings.warn(f"Found negative wavenumbers: '{neg_nus}'!")

    geom_redund = geom.copy_all(
        coord_type="redund",
        coord_kwargs={
            "bonds_only": True,
        },
    )
    B = geom_redund.internal.B
    typed_prims = geom_redund.internal.typed_prims

    return stretch_constants(hessian, L, typed_prims, B)
