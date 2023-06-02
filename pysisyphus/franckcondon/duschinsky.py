# [1] https://pubs.acs.org/doi/pdf/10.1021/jp004230b
#     Ab Initio Computation of the Duschinsky Mixing of Vibrations and
#     Nonlinear Effects
#     Sando, Spears, 2001
# [2] https://pubs.acs.org/doi/pdf/10.1021/jp004229c
#     Large Electron Transfer Rate Effects from the Duschinsky Mixing of Vibrations
#     Sando, Spears, 2001
#     Currently, nothing from [2] is implemented here
# [3] https://pubs.acs.org/doi/pdf/10.1021/acs.jctc.7b00325
#     Franck−Condon Models for Simulating the Band Shape of
#     Electronic Absorption Spectra
#     Li, Truhlar, 2017


from dataclasses import dataclass
from functools import singledispatch
from enum import Enum
import itertools as it

import numpy as np
import scipy as sp

from pysisyphus.constants import AMU2AU, AMU2KG, BOHR2M, C, HBAR
from pysisyphus.Geometry import Geometry


AMU2KG_SQRT = np.sqrt(AMU2KG)


@dataclass
class AxisSwitch:
    T0: np.ndarray
    B0: np.ndarray
    masses: np.ndarray
    coords3d_init: np.ndarray
    coords3d_final: np.ndarray
    coords3d_init_rot: np.ndarray


def get_axis_switch(coords3d_init, coords3d_final, masses, planar=False):
    assert coords3d_init.shape == coords3d_final.shape
    natoms, _ = coords3d_init.shape

    coords_init = coords3d_init.flatten()
    coords_final = coords3d_final.flatten()

    # TODO: handle planarity; set C_zz unity
    if planar:
        raise Exception("Implement handling of planar coordinates!")
    # Eq. (11) in [1]
    C = np.einsum("i,ia,ib->ab", masses, coords3d_final, coords3d_init)
    C_inv = np.linalg.pinv(C)
    eigvals, R = np.linalg.eigh(C.T @ C)
    eigval_mat = np.diag(np.sqrt(eigvals))
    ones = (-1, 1)

    def make_B0(T0):
        return sp.linalg.block_diag(*[T0 for _ in range(natoms)])

    best_diff_norm = None
    for ones_diag in it.product(ones, ones, ones):
        # Either -1 or 1 on diagonal; 2³ = 8 different matrices possible.
        lambda_mat = np.diag(ones_diag)
        # Eq. (12) in [1]
        T0 = R @ lambda_mat @ eigval_mat @ R.T @ C_inv
        det = np.linalg.det(T0)
        # det(T0) must be +1
        if abs(det - 1) > 1e-14:
            continue
        B0 = make_B0(T0)
        # Be careful to use B.T, not B. Here, we basically apply eq. (10) in [1].
        # In the paper (T⁰)⁻¹ is used. As T⁰ is orthogonal, the inverse is its
        # transpose.
        coords_init_rot = B0.T @ coords_init
        diff = coords_init_rot - coords_final
        diff_norm = np.linalg.norm(diff)
        if (best_diff_norm is None) or (diff_norm < best_diff_norm):
            best_diff_norm = diff_norm
            best_axis_switch = AxisSwitch(
                T0,
                B0,
                masses,
                coords3d_init,
                coords3d_final,
                coords_init_rot.reshape(-1, 3),
            )
    return best_axis_switch


DuschinskyRef = Enum("DuschinskyRef", ("INITIAL", "FINAL"))


@dataclass
class DuschinskyResult:
    J: np.ndarray
    K: np.ndarray
    ref: DuschinskyRef

    def to_K_unitless(self, wavenums):
        # Wavenumbers in m⁻¹
        wavenums_m = wavenums * 1e2
        # First term results in kg Bohr². We have to transform it to sqrt(AMU) * Bohr²
        conv = (
            1 / np.sqrt(HBAR / (2 * np.pi * wavenums_m * C * BOHR2M**2)) * AMU2KG_SQRT
        )
        K_unitless = self.K * conv
        return K_unitless


@singledispatch
def duschinsky(
    L_init: np.ndarray,
    coords3d_init: np.ndarray,
    L_final: np.ndarray,
    coords3d_final: np.ndarray,
    masses: np.ndarray,
    reference: DuschinskyRef = DuschinskyRef.INITIAL,
    with_axis_switch=True,
) -> DuschinskyResult:
    """Duschinsky matrix & displacment matrix according to [1].

    Q' = JQ + K
    J = L'^T L  (Duschinsky matrix J)
    K = L'^T m^(1/2) * (x^0 - x^0') (Displacement vector)
    """

    # Shape (3N, nmodes)
    assert L_init.shape == L_final.shape
    _3natoms, _ = L_init.shape
    coords_init = coords3d_init.flatten()
    coords_final = coords3d_final.flatten()
    assert coords_init.size == _3natoms
    assert coords_final.size == _3natoms
    assert masses.size == _3natoms // 3

    # Displacements in terms of final normal coordinates
    if reference == DuschinskyRef.FINAL:
        L = L_final
        coords_init, coords_final = coords_final, coords_init
    # Displacements in terms of initial normal coordinates
    elif reference == DuschinskyRef.INITIAL:
        L = L_init
    else:
        raise Exception("Invalid reference!")

    if with_axis_switch:
        # def get_axis_switching_B(coords3d_init, coords3d_final, masses, planar=False):
        axis_switch = get_axis_switch(coords3d_init, coords3d_final, masses)
        B0 = axis_switch.B0
    else:
        B0 = np.eye(_3natoms)

    # Duschinsky matrix
    J = L_final.T @ B0.T @ L_init
    M = np.sqrt(np.repeat(masses, 3))
    coords_diff = B0.T @ coords_final - coords_init
    K = L.T * M @ coords_diff

    dusch_res = DuschinskyResult(J=J, K=K, ref=reference)
    return dusch_res


@duschinsky.register
def _(geom_init: Geometry, geom_final: Geometry, **kwargs):
    # nus, eigvals, mw_cart_displs, cart_displs
    _, _, L_init, _ = geom_init.get_normal_modes()
    coords3d_init = geom_init.coords3d
    _, _, L_final, _ = geom_final.get_normal_modes()
    coords3d_final = geom_final.coords3d
    masses = geom_init.masses
    return duschinsky(L_init, coords3d_init, L_final, coords3d_final, masses, **kwargs)


def unitless_displs_from_eigensystem(
    mw_gradient: np.ndarray, eigenvalues: np.ndarray, eigenvectors: np.ndarray
):
    """Unitless displacements from gradient and normal mode eigensystem.

    Useful to determine Huang-Rhys factors and/or displacements in situations
    where no Hessian is avaialable or its calculation is impossible, e.g., in
    excited state calculations.

    Parameters
    ----------
    mw_gradient
        Mass-weighted gradient in atomic units Eh/(a0 sqrt(amu)).
    eigenvalues
        Eigenvalues of the projected, mass-weighted Hessian in Eh/(a0² amu).
    eigenvectors
        Eigenvectors of the projected, mass-weighted Hessian. Unitless.

    Returns
    -------
    Unitless displacements along the normal modes.
    """
    A = eigenvectors
    inv_eigvals = 1 / eigenvalues
    # Eq. (16) of [3]
    q = -A @ np.diag(inv_eigvals) @ A.T @ mw_gradient
    # Eq. (19) of [3]
    q_tilde = A.T @ q  # Unit of sqrt(amu) * Bohr

    # (ME/AMU)**(1/4)
    conv = (1 / AMU2AU) ** 0.25
    qsqrt_eigvals = eigenvalues**0.25 / conv
    # Eq. (20) of [3]. If this seems arcane to you, write down the units,
    # then you'll see it.
    delta = np.diag(qsqrt_eigvals) @ q_tilde
    return delta
