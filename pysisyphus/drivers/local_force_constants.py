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


import numpy as np


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
    """
    Parameters
    ----------
    hessian
        Cartesian Hessian matrix.
    B
        Wilson-B-matrix.
    L
        Renormalized Cartesian displacment vectors. Obtained from unweighting
        normal modes and renormalizing.

    Returns
    -------
    force_constants
        Local force constants in atomic units (mass/time**2).
    """

    # Eq. (10) in [1]
    D = B @ L
    # Eq. (6) in [1]
    K = L.T @ hessian @ L
    K_inv = np.diag(1/np.diag(K))
    force_constants = np.zeros(B.shape[0])
    for i, d_row in enumerate(D):
        force_constants[i] = 1 / (d_row @ K_inv @ d_row.T)
    return force_constants
