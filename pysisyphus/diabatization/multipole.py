# [1] https://doi.org/10.1063/1.4894472
#     Diabatization based on the dipole and quadrupole: The DQ method
#     Hoyer, Xu, Ma, Gagliardi, Truhlar, 2014
# [2] https://doi.org/10.1063/1.4948728
#     The DQ and DQΦ electronic structure diabatization methods:
#     Validation for general applications
#     Hoyer, Parker, Galgiardi, Truhlar, 2016
# [3] https://doi.org/10.1063/1.3042233
#     Constructing diabatic states from adiabatic states: Extending
#     generalized Mulliken–Hush to multiple charge centers with Boys localization
#     Subotnik, Yeganeh, Cave, Ratner, 2008
# [4]  https://doi.org/10.1007/BF01113521
#      Efficient use of Jacobi rotations for orbital optimization and localization
#      Raffenetti, Ruedenberg, Janssen, Schaefer, 1993
# [5]  https://doi.org/10.1103/RevModPhys.35.457
#      Localized Atomic and Molecular Orbitals
#      Edmiston, Ruedenberg, 1963
# [6] https://doi.org/10.1063/1.3148777
#      The initial and final states of electron and energy transfer processes:
#      Diabatization as motivated by system-solvent interactions
#      Subotnik, Cave, Steele, Shenvi, 2009

from typing import Optional

from jinja2 import Template
import numpy as np

from pysisyphus.diabatization.helpers import get_random_U
from pysisyphus.diabatization.results import (
    DiabatizationResult,
    dia_result_from_jac_result,
)
from pysisyphus.helpers_pure import highlight_text
from pysisyphus.wavefunction import logger
from pysisyphus.wavefunction.localization import (
    JacobiSweepResult,
    jacobi_sweeps,
    get_fb_ab_func,
    get_fb_cost_func,
)


DQTemplate = Template(
    """{{ name }}

Dipole moments
--------------
{{ dip_moms }}
{%- if quad_moms is not none%}
Trace of quadrupole tensor
--------------------------
{{ quad_moms }}
α = {{ "%.2f"|format(alpha) }}
{% endif %}
{%- if epots is not none %}
Electronic component of electrostatic potential
-----------------------------------------------
{{ epots}}
β = {{ "%.2f"|format(beta) }}
{% endif %}

"""
)


def dq_jacobi_sweeps(
    dip_moms: np.ndarray,
    quad_moms: Optional[np.ndarray] = None,
    epots: Optional[np.ndarray] = None,
    alpha: Optional[float] = 10.0,
    beta: Optional[float] = 1.0,
    random: bool = False,
) -> JacobiSweepResult:
    """Rotation matrix from DQ-diabatization as outlined in [1], [2] and [3].

    When no quadrupole moments are given, the DQ-diabatization reduces to a simple
    Boys-diabatization, as outlined by Subotnik et al in [3]. In this case we just
    zeros for the quadrupole moments. As the overall size of the matrices is small,
    the additional FLOPs dont hurt and the code can be kept simpler.

    We only use the trace of the quadrupole moment matrix. There are three
    dipole moment components, but only one trace of the quadrupole moment
    matrix.

    TODO: In principle this function can easily be extended to support an arbitrary
    number of properties, each with its own scaling factor. When this is implemented,
    we could drop the separate Boys-localiatzion function in wavefunction.localization.

    Parameters
    ----------
    dip_moms
        Dipole moment matrix of the adiabatic states.
        Shape (3, nstates, nstates).
    quad_moms
        Optional matrix containing the trace of the primitive quadrupole
        moments. Shape (1, nstates, nstates). Optional.
    epots
        Electronic part of the electrostatic potential.
    alpha
        Scaling factor for quadrupole moment contribution.
    beta
        Scaling factor for electrostatic potential contribution.
    random
        Boolean that controls if we start from the original adiabatic states (rotation
        matrix U is the identity matrix) or if we start from randomly mixed states. In
        high symmetry systems it may be beneficial to start from a random state.

    Returns
    -------
    JacobiSweepResult
        Diabatization result.
    """

    _, nstates, _ = dip_moms.shape

    if random:
        U = get_random_U(nstates)
    else:
        U = np.eye(nstates)

    assert dip_moms.shape == (3, *U.shape)
    assert U.ndim == 2
    assert U.shape[0] == U.shape[1]
    expected_quad_mom_shape = (1, *dip_moms.shape[1:])

    prefixes = [
        char_
        for property_, char_ in ((dip_moms, "D"), (quad_moms, "Q"), (epots, "Φ"))
        if property_ is not None
    ]
    name = "".join(prefixes) + "-diabatization"
    if no_quad_moms := quad_moms is None:
        quad_moms = np.zeros(expected_quad_mom_shape)
    if no_epots := epots is None:
        epots = np.zeros(expected_quad_mom_shape)
    quad_moms = quad_moms.reshape(*expected_quad_mom_shape)
    epots = epots.reshape(*expected_quad_mom_shape)

    dip_ab_func = get_fb_ab_func(dip_moms)
    quad_ab_func = get_fb_ab_func(quad_moms)
    epot_ab_func = get_fb_ab_func(epots)

    def ab_func(j, k, U):
        dip_A, dip_B = dip_ab_func(j, k, U)
        quad_A, quad_B = quad_ab_func(j, k, U)
        epot_A, epot_B = epot_ab_func(j, k, U)
        A = dip_A + alpha * quad_A + beta * epot_A
        B = dip_B + alpha * quad_B + beta * epot_B
        return A, B

    dip_cost_func = get_fb_cost_func(dip_moms)
    quad_cost_func = get_fb_cost_func(quad_moms)
    epot_cost_func = get_fb_cost_func(epots)

    def cost_func(U):
        dip_P = dip_cost_func(U)
        quad_P = quad_cost_func(U)
        epot_P = epot_cost_func(U)
        P = dip_P + alpha * quad_P + beta * epot_P
        return P

    msg = DQTemplate.render(
        name=highlight_text(name),
        dip_moms=dip_moms,
        quad_moms=quad_moms if not no_quad_moms else None,
        alpha=alpha,
        epots=epots if not no_epots else None,
        beta=beta,
    )
    logger.info(msg)

    return jacobi_sweeps(U, cost_func, ab_func)


def dq_diabatization(
    adia_ens: np.ndarray,
    dip_moms: np.ndarray,
    quad_moms: Optional[np.ndarray] = None,
    epots: Optional[np.ndarray] = None,
    **kwargs,
) -> DiabatizationResult:
    """Property-based diabatization using multipole moments.

    Similar to Foster-Boys-localization, but w/ electronic states.

    Parameters
    ----------
    adia_ens
        1d array of shape (nstates, ) containing adiabatic excitation energies
        in atomic units/Hartree.
    dip_moms
        3d array of shape (3, nstates, nstates) containing (transition) dipole moments
        (x, y, z).
    quad_moms
        Optional 3d array of shape (3, nstates, nstates) containing the diagonal of
        the (transition) quadrupole moment tensor (xx, yy, zz)
    epots
        Optional 3d array of shape (3, nstates, nstates) containing electronic part
        of electrostatic potential.

    Returns
    -------
    dia_result
        Result of diabatization containing various quantities, e.g., the adiabatic-
        diabatic-transformation matrix.
    """
    jac_res = dq_jacobi_sweeps(dip_moms, quad_moms=quad_moms, epots=epots, **kwargs)
    kind = "dq"
    return dia_result_from_jac_result(
        kind, adia_ens, jac_res, dip_moms=dip_moms, quad_moms=quad_moms, epots=epots
    )
