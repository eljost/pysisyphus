from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

from pysisyphus.constants import AU2EV
from pysisyphus.wavefunction import logger
from pysisyphus.wavefunction.localization import (
    JacobiSweepResult,
    jacobi_sweeps,
    get_fb_ab_func,
    get_fb_cost_func,
)


def dq_jacobi_sweeps(
    dip_moms: NDArray[float],
    quad_moms: Optional[NDArray[float]] = None,
    alpha: Optional[float] = 1.0,
) -> JacobiSweepResult:
    """Rotation matrix from DQ-diabatization as outlined in [5].

    When no quadrupole moments are given, the DQ-diabatization reduces to a simple
    Boys-diabatization, as outlined by Subotnik et al in [3]. In this case we just
    zeros for the quadrupole moments. As the overall size of the matrices is small,
    the additional FLOPs dont hurt and the code can be kept simpler.

    We only use the trace of the quadrupole moment matrix. There are three
    dipole moment components, but only one trace of the quadrupole moment
    matrix.

    Parameters
    ----------
    dip_moms
        Dipole moment matrix of the adiabatic states.
        Shape (3, nstates, nstates).
    quad_moms
        Optional matrix containing the trace of the primitive quadrupole
        moments. Shape (1, nstates, nstates). Optional.
    alpha
        Factor controlling the quadrupole moment contribution.


    Returns
    -------
    JacobiSweepResult
        Diabatization result.
    """

    _, nstates, _ = dip_moms.shape
    U = np.eye(nstates)

    assert dip_moms.shape == (3, *U.shape)
    assert U.ndim == 2
    assert U.shape[0] == U.shape[1]
    expected_quad_mom_shape = (1, *dip_moms.shape[1:])
    if quad_moms is None:
        name = "Boys-diabatization"
        quad_moms = np.zeros(expected_quad_mom_shape)
    else:
        name = "DQ-diabatization"
    quad_moms = quad_moms.reshape(*expected_quad_mom_shape)

    dip_ab_func = get_fb_ab_func(dip_moms)
    quad_ab_func = get_fb_ab_func(quad_moms)

    def ab_func(j, k, U):
        dip_A, dip_B = dip_ab_func(j, k, U)
        quad_A, quad_B = quad_ab_func(j, k, U)
        A = dip_A + alpha * quad_A
        B = dip_B + alpha * quad_B
        return A, B

    dip_cost_func = get_fb_cost_func(dip_moms)
    quad_cost_func = get_fb_cost_func(quad_moms)

    def cost_func(U):
        dip_P = dip_cost_func(U)
        quad_P = quad_cost_func(U)
        P = dip_P + alpha * quad_P
        return P

    logger.info(name)
    return jacobi_sweeps(U, cost_func, ab_func)


def dq_diabatization_from_npz(npz_fn, use_states=None):
    data = np.load(npz_fn)
    energies = data["energies"]
    dipoles = data["dipoles"]
    quadrupoles = data["quadrupoles"]

    energies -= energies.min()
    energies *= AU2EV

    if use_states is not None:
        energies = energies[:, use_states]

    adia_energies = np.zeros_like(energies)
    dia_energies = np.zeros_like(energies)
    for i, (ens, dpm, qpm) in enumerate(zip(energies, dipoles, quadrupoles)):
        tr_qpm = qpm[0, 0] + qpm[1, 1] + qpm[2, 2]
        dia_res = dq_jacobi_sweeps(dpm, tr_qpm, alpha=10.0)
        U = dia_res.C
        print(f"Rotation matrix\n{U}\n")

        adia_energies[i] = ens
        adia_mat = np.diag(ens)
        dia_mat = U.T @ adia_mat @ U
        dia_ens = np.diag(dia_mat)
        dia_energies[i] = dia_ens
    return adia_energies, dia_energies


def plot_adia_dia(adia_energies, dia_energies, show=True):
    fig, ax = plt.subplots()
    for i, state in enumerate(adia_energies.T):
        ax.plot(state, "o--", label=f"$V_{i}$")
    for i, state in enumerate(dia_energies.T):
        ax.plot(state, "x-", label=f"$U_{i}$")
    ax.legend()
    ax.set_xlabel("Step")
    ax.set_ylabel("$\Delta E$ / eV")
    if show:
        plt.show()
    return fig, ax
