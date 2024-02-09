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

from dataclasses import dataclass
from typing import Optional

from jinja2 import Template
import numpy as np
from numpy.typing import NDArray
from scipy.stats import special_ortho_group

from pysisyphus.helpers_pure import highlight_text, to_subscript_num
from pysisyphus.wavefunction import logger
from pysisyphus.wavefunction.localization import (
    JacobiSweepResult,
    jacobi_sweeps,
    get_fb_ab_func,
    get_fb_cost_func,
)


def get_random_U(N):
    """Get random rotation matrix."""
    return special_ortho_group.rvs(N)


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
    dip_moms: NDArray[float],
    quad_moms: Optional[NDArray[float]] = None,
    epots: Optional[NDArray[float]] = None,
    alpha: Optional[float] = 10.0,
    beta: Optional[float] = 1.0,
    random: bool = False,
    jc_kwargs: Optional[dict] = None,
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

    if jc_kwargs is None:
        jc_kwargs = {}

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

    return jacobi_sweeps(U, cost_func, ab_func, **jc_kwargs)


ERTemplate = Template(
    """{{ name }}

Coulomb-Tensor R
--------------
{{ R }}

"""
)


def edmiston_ruedenberg_jacobi_sweeps(
    R: np.ndarray, random: bool = False, jc_kwargs: Optional[dict] = None
) -> JacobiSweepResult:
    """Edmiston-Ruedenberg localization using Jacobi sweeps.

    By any means, this is not intended to be used to localize (many) orbitals,
    but to diabatize adiabatic electronic states. The tensor R contains the
    expectation values of the 2-electron Coulomb operator for the different
    adiabatic states and the transition densities between them.

    The rotated tensor R_rot is a quartic function of the 2d rotation matrix U and
    the current implementation using einsum scales as N⁵, with N being the number
    of states. As the this function is intended to be used for diabatization, N
    should be a small number, so the N⁵ shouldn't matter.

    Parameters
    ----------
    R
        4d-tensor of shape (nstates, nstates, nstates, nstates), containing
        2-electron integrals contracted with the relaxed
        densities of the different states and the transition densities between
        them. See eqs. (A8) and (A10) in [6] on how to calculate it.
    random
        Boolean that controls if we start from the original adiabatic states (rotation
        matrix U is the identity matrix) or if we start from randomly mixed states. In
        high symmetry systems it may be beneficial to start from a random state.
    jc_kwargs
        Optinal kwargs for the Jacobi-Sweeps function.

    Returns
    -------
    JacobiSweepResult
        Dataclass containing all relevant parameters & results of the Jacobi-Sweeps
        run.
    """
    if jc_kwargs is None:
        jc_kwargs = {}

    nstates = len(R)
    assert R.shape == (nstates, nstates, nstates, nstates)

    if random:
        U = get_random_U(nstates)
    else:
        U = np.eye(nstates)

    def ab_func(j, k, U):
        # See eq. (2.2) in [4]
        R_rot = np.einsum("jJ,kK,lL,mM,JKLM->jklm", U, U, U, U, R, optimize="greedy")
        # The attentive reader may notice, that the implementation of A below is not
        # eq. (17) from [5] but it is the negative of it. Depending on how the angle α
        # is calculated (three variants are given in eq. (19) in [5]), different arguments
        # are used for the trigonometric functions. Our Jacobi-sweeps implementation
        # uses arccos, to recover the angle from, so we have to return -A, instead of A.
        # The sign of B doesn't matter, as it only appears squared in the denominator of
        # eq. (19) in [5].
        A = -R_rot[j, k, j, k] + 0.25 * (
            R_rot[j, j, j, j] - 2 * R_rot[j, j, k, k] + R_rot[k, k, k, k]
        )
        B = R_rot[j, j, j, k] - R_rot[k, k, j, k]
        return A, B

    def cost_func(U):
        # Eq. (A9) in [6]
        return np.einsum("JI,KI,LI,MI,JKLM->", U, U, U, U, R)

    msg = ERTemplate.render(
        name=highlight_text("Edmiston-Ruedenberg-diabatization"),
        R=R,
    )
    logger.info(msg)

    return jacobi_sweeps(U, cost_func, ab_func, **jc_kwargs)


DiaResultTemplate = Template(
    """
All energies are given in {{ unit }}.

Adiabatic energy matrix V
-------------------------
{{ adia_mat }}

Rotation matrix U
-----------------
{{ U }}
det(U)={{ "%.4f"|format(det) }}

Diabatic energy matrix D = UᵀVU
-------------------------------
{{ dia_mat }}

Diabatic states Ξᵢ sorted by energy
--------------------------------
{%- for ind, dia_state, dia_en in dia_states_sorted %}
{{ ind }}: {{ dia_state }}, {{ "%.4f"|format(dia_en) }} {{ unit }}
{%- endfor %}

Composition of diabatic states Ξᵢ
---------------------------------
{%- for dia_comp in dia_compositions %}
{{ dia_comp }}
{%- endfor %}

Weights U²
----------
{{ U2 }}

Absolute diabatic couplings
---------------------------
{%- for key, coupling in couplings %}
|{{ key }}| = {{ "%.4f"|format(coupling) }} {{ unit }}
{%- endfor %}
"""
)


@dataclass
class DiabatizationResult:
    U: np.ndarray
    adia_ens: np.ndarray
    # U and adia_ens should be enough ... the rest should be optional.
    # dia_ens can be calculated by rotationg the adiabatic energy matrix,
    # which is just a diagonal matrix containing the adiabatic energies.
    dia_ens: np.ndarray
    adia_mat: np.ndarray
    dia_mat: np.ndarray
    is_converged: bool
    cur_cycle: int
    P: float
    # TODO: store rotated/mixed properties instead of original ones?
    dip_moms: Optional[np.ndarray] = None
    quad_moms: Optional[np.ndarray] = None
    epots: Optional[np.ndarray] = None
    R: Optional[np.ndarray] = None

    @property
    def nstates(self):
        return len(self.adia_ens)

    def __post_init__(self):
        assert len(self.adia_ens) == len(self.dia_ens)
        quad_shape = (self.nstates, self.nstates)
        assert self.U.shape == self.adia_mat.shape == self.dia_mat.shape == quad_shape

    @property
    def couplings(self):
        nstates = self.nstates
        abs_dia_mat = np.abs(self.dia_mat)
        couplings = dict()
        for from_ in range(nstates):
            for to_ in range(from_ + 1, nstates):
                key = (from_, to_)
                couplings[key] = couplings[key[::-1]] = abs_dia_mat[from_, to_]
        return couplings

    def render_report(self, adia_labels=None, unit="eV"):
        U = self.U
        nstates = self.nstates
        U2 = U**2
        det = np.linalg.det(U)

        subscripts = [to_subscript_num(i) for i in range(nstates)]
        adia_wfs = [f"Φ{subscripts[i]}" for i in range(nstates)]
        dia_wfs = [f"Ξ{subscripts[i]}" for i in range(nstates)]

        # Diabatic states, sorted by energy
        dia_ens = self.dia_ens
        sorted_inds = np.argsort(dia_ens)
        dia_states_sorted = list()
        for i, sort_ind in enumerate(sorted_inds):
            dia_states_sorted.append((i, dia_wfs[sort_ind], dia_ens[sort_ind]))

        # Composition of diabatic states
        dia_compositions = list()
        if adia_labels:
            adia_labels = [f"({al})" for al in adia_labels]
        else:
            adia_labels = ["" for _ in range(nstates)]
        for i in range(nstates):
            coeffs = U[:, i]
            signs = ["+" if c > 0 else "-" for c in coeffs]
            sum_ = " ".join(
                [
                    f"{sign} {c:>6.4f}·{awf}{al}"
                    for sign, c, awf, al in zip(
                        signs, np.abs(coeffs), adia_wfs, adia_labels
                    )
                ]
            )
            dia_compositions.append(f"{dia_wfs[i]} = {sum_}")

        # Couplings
        couplings = self.couplings
        flat_couplings = list()
        for from_ in range(nstates):
            for to_ in range(from_ + 1, nstates):
                key = f"D{subscripts[from_]}{subscripts[to_]}"
                flat_couplings.append((key, couplings[(from_, to_)]))

        rendered = DiaResultTemplate.render(
            unit=unit,
            adia_mat=self.adia_mat,
            U=U,
            det=det,
            dia_mat=self.dia_mat,
            dia_states_sorted=dia_states_sorted,
            dia_compositions=dia_compositions,
            U2=U2,
            couplings=flat_couplings,
        )
        return rendered


def dq_diabatization(adia_ens, dip_moms, quad_moms=None, epots=None, **kwargs):
    jac_res = dq_jacobi_sweeps(dip_moms, quad_moms=quad_moms, epots=epots, **kwargs)
    U = jac_res.C
    adia_mat = np.diag(adia_ens)
    dia_mat = U.T @ adia_mat @ U
    dia_ens = np.diag(dia_mat)
    dia_result = DiabatizationResult(
        U=U,
        adia_ens=adia_ens,
        dia_ens=dia_ens,
        adia_mat=adia_mat,
        dia_mat=dia_mat,
        is_converged=jac_res.is_converged,
        cur_cycle=jac_res.cur_cycle,
        P=jac_res.P,
        dip_moms=dip_moms,
        quad_moms=quad_moms,
        epots=epots,
    )
    return dia_result


def edmiston_ruedenberg_diabatization(adia_ens, R, **kwargs):
    jac_res = edmiston_ruedenberg_jacobi_sweeps(R, **kwargs)
    U = jac_res.C
    adia_mat = np.diag(adia_ens)
    dia_mat = U.T @ adia_mat @ U
    dia_ens = np.diag(dia_mat)
    dia_result = DiabatizationResult(
        U=U,
        adia_ens=adia_ens,
        dia_ens=dia_ens,
        adia_mat=adia_mat,
        dia_mat=dia_mat,
        is_converged=jac_res.is_converged,
        cur_cycle=jac_res.cur_cycle,
        P=jac_res.P,
    )
    return dia_result
