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
import itertools as it
from typing import Optional
import warnings

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
    edmiston_ruedenberg,
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


ERTemplate = Template(
    """{{ name }}

Coulomb-Tensor R
----------------
{{ R }}

"""
)


def edmiston_ruedenberg_jacobi_sweeps(
    R: np.ndarray,
    random: bool = False,
    U=None,
    max_cycles: int = None,
    dP_thresh: float = 1e-8,
) -> JacobiSweepResult:
    assert R.ndim == 4
    nstates = R.shape[0]
    assert all([dim == nstates for dim in R.shape])

    if max_cycles is None:
        max_cycles = nstates * 100

    if U is None:
        if random:
            U = get_random_U(nstates)
        else:
            U = np.eye(nstates)

    msg = ERTemplate.render(
        name=highlight_text("Edmiston-Ruedenberg-diabatization"),
        R=R,
    )
    logger.info(msg)

    D_prev = np.nan
    # Macro-cycles
    for i in range(max_cycles):
        # Calculate rotated Coulomb-tensor
        #
        # The rotated Coulomb tensor is a quartic function of the rotation matrix 'U'.
        # The exact contraction of the Coulomb tensor 'R' with 'U' depends on whether rows
        # or columns of 'U' are rotated/mixed.
        # When rows of 'U' are mixed/rotated we have to use
        #   R_rot = np.einsum("jJ,kK,lL,mM,JKLM->jklm", U, U, U, U, R, optimize="greedy") .
        # When columns of 'U' are mixed/rotated we have to use
        R_rot = np.einsum("Jj,Kk,Ll,Mm,JKLM->jklm", U, U, U, U, R, optimize="greedy")
        #
        # Even though this function is intended to be used for diabatization where 'U' will
        # be quiet small we chose to rotate columns of 'U'. This way, this function can also
        # be used for MO-localization, where the different MOs will be organized in columns.
        # For the diabatization of electronic states it does not matter if we rotate rows or
        # columns, as long as the rotation and the calculation of R_rot are done in a consistent
        # way.
        # The calculation of R_rot scales with nstates**5, but this will be no problem when
        # nstates remains small.

        # Edmiston-Ruedenberg cost-function
        D = np.einsum("IIII->", R_rot)
        # Cost-function change compared to previous macro-cycle
        dD = D - D_prev
        print(f"Macro cycle {i: >3d}: {D=: >14.6f}, {dD=: >+12.6e}")
        if converged := (abs(dD) <= dP_thresh):
            print("Converged!")
            break
        D_prev = D

        # Micro cycles
        #
        # Loop over all possible state pairs and determine the pair that promises the
        # greatest cost-function increase.
        term_max = None
        j_max = None
        k_max = None
        A_max = None
        B_max = None
        for j, k in it.combinations(range(nstates), 2):
            R1212 = R_rot[j, k, j, k]
            R1111 = R_rot[j, j, j, j]
            R1122 = R_rot[j, j, k, k]
            R2222 = R_rot[k, k, k, k]
            R1112 = R_rot[j, j, j, k]
            R1222 = R_rot[j, k, k, k]

            # Eq. (19) in [5
            A = R1212 - 0.25 * (R1111 - 2.0 * R1122 + R2222)
            B = R1112 - R1222
            # RHS of eq. (26) in [5]
            term = A + (A**2 + B**2) ** 0.5
            # print(
            # f"\t{j=: >3d}, {k=: >3d}, {A=: >12.6f}, {B=: >12.6f}, {term=: >12.6f}"
            # )
            # Update pair that promises the greatest cost-function increase
            if (term_max is None) or term > term_max:
                j_max = j
                k_max = k
                term_max = term
                A_max = A
                B_max = B
        # End loop over all micro-cycles

        # When B is very small and A is negative, then A + (A**2 + B**2)**0.5 will be 0.0
        # and term_max will be 0.0 too. In this case no rotation will improve the cost-function
        # to an appreciable degree and we skip the rotation. This will result in convergence
        # in the next macro-cycle, as the cost-function will stay constant.
        if abs(term_max) <= 1e-12:
            continue

        j = j_max
        k = k_max
        A = A_max
        B = B_max
        if term_max is None:
            print("A**2 + B**2 is very small. Indicating convergence.")
            break
        # Denominator sqrt-term of eq. (19) in [5]
        sqrt_term = (A**2 + B**2) ** 0.5
        # Determine rotation angle according to right column on p. 460 of [5]
        # Calcuate cosine and sinus terms according to eq. (19) in [5]
        cos4a = -A / sqrt_term
        sin4a = B / sqrt_term
        # Two pairs (x1, y1) and (x2, y2)
        #
        # First pair
        x21 = 0.5 * (1.0 + (1.0 - 0.5 * (1 - cos4a)) ** 0.5)
        #                ^--- Note this sign
        x1 = x21**0.5
        y1 = (1.0 - x21) ** 0.5
        # Second pair
        x22 = 0.5 * (1.0 - (1.0 - 0.5 * (1 - cos4a)) ** 0.5)
        #                ^--- Note this sign
        x2 = x22**0.5
        y2 = (1.0 - x22) ** 0.5
        # Check which pair fulfills 4 * x * y * (x**2 - y**2) = sin(4a)
        for x, y in ((x1, y1), (x2, y2)):
            diff = 4 * x * y * (x**2 - y**2) - sin4a
            # Interestingly, diff can become quite 'big' sometimes (> 1e-8)
            if abs(diff) <= 5e-8:
                break
        else:
            raise Exception(
                f"Determination of rotation angle failed ({diff=: >12.6e})!"
            )
        # x of the correct pair corresponds to the cosine of the rotation angle
        rad = np.arccos(x)
        # deg = np.rad2deg(rad)
        # print(f"\tRotation of {j=} and {k=} by {rad=: 10.6e} rad ({deg: >10.6f}°).")
        # Inplace rotation of matrix columns with indices 'j' and 'k' by 'rad' radians.
        # Eq. (15) in [5]
        # NOTE: whether we rotate/mix columns or rows of U determines how we have to rotate
        # the Coulomb-tensor at the beginning of a macro-cycle.
        cos = np.cos(rad)
        sin = np.sin(rad)
        j_rot = cos * U[:, j] + sin * U[:, k]
        k_rot = -sin * U[:, j] + cos * U[:, k]
        U[:, j] = j_rot
        U[:, k] = k_rot
    else:
        raise Exception(f"ER-localization did not converge after {max_cycles} cycles!")
    # End loop over macro-cycles

    result = JacobiSweepResult(
        is_converged=converged,
        cur_cycle=i,
        C=U,
        P=D,
    )
    return result


DiaResultTemplate = Template(
    """
########################
# DIABATIZATION REPORT #
########################

All energies are given in {{ unit }}.

Every column of the rotation matrix U describes the composition of a
diabatic state in terms of (possibly various) adiabatic states.

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
-----------------------------------
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

Unique absolute diabatic couplings
----------------------------------
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
    # Dipole moments
    dip_moms: Optional[np.ndarray] = None
    # Quadrupole moments
    quad_moms: Optional[np.ndarray] = None
    # Electrostatic potentials
    epots: Optional[np.ndarray] = None
    # 4d Coulomb tensor  of shape (nstates, nstates, nstates, nstates)
    R_tensor: Optional[np.ndarray] = None  # Coulomb tensor
    # 3d density fitting tensor of shape (naux, nstates, nstates)
    L_tensor: Optional[np.ndarray] = None
    # TODO: add adiabatic labels and use them in the report

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


def dia_result_from_jac_result(
    adia_ens: np.ndarray, jac_res: JacobiSweepResult, **property_tensors
) -> DiabatizationResult:
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
        **property_tensors,
    )
    return dia_result


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
    return dia_result_from_jac_result(
        adia_ens, jac_res, dip_moms=dip_moms, quad_moms=quad_moms, epots=epots
    )


def edmiston_ruedenberg_diabatization_jacobi(
    adia_ens: np.ndarray, R_tensor: np.ndarray, **kwargs
) -> DiabatizationResult:
    """Property-based diabatization using the Coulomb tensor and Jacobi sweeps.

    Similar to Edmiston-Ruedenberg-localization, but w/ electronic states.

    Parameters
    ----------
    adia_ens
        1d array of shape (nstates, ) containing adiabatic excitation energies
        in atomic units/Hartree.
    R_tensor
        4d array of shape (nstates, nstates, nstates, nstates). Coulomb tensor
        in the basis of the adiabatic electronic states.

    Returns
    -------
    dia_result
        Result of diabatization containing various quantities, e.g., the adiabatic-
        diabatic-transformation matrix.
    """
    jac_res = edmiston_ruedenberg_jacobi_sweeps(R_tensor, **kwargs)
    return dia_result_from_jac_result(
        adia_ens,
        jac_res,
        R_tensor=R_tensor,
    )


# Shortcut
edmiston_ruedenberg_diabatization = edmiston_ruedenberg_diabatization_jacobi


def edmiston_ruedenberg_diabatization_df(
    adia_ens: np.ndarray,
    df_tensor: np.ndarray,
    overlap_matrix: Optional[np.ndarray] = None,
    nruns: int = 25,
    max_cycles: int = 10_000,
) -> DiabatizationResult:
    """Property-based diabatization using the Coulomb tensor and density fitting.

    Similar to Edmiston-Ruedenberg-localization, but w/ electronic states.

    TODO: factor out the call to edmiston_ruedenberg() and put it in a separate function.

    Parameters
    ----------
    adia_ens
        1d array of shape (nstates, ) containing adiabatic excitation energies
        in atomic units/Hartree.
    df_tensor
        3d array of shape (naux, nstates, nstates). Coulomb tensor
        in the basis of the adiabatic electronic states.
    overlap_matrix
        Optional 2d array of shape (nstates, nstates). If not given, orthogonal
        electronic states are assumed, which should be the case in typical TDDFT
        calculations.
    nruns
        Number of initial conditions that are tested/number of macro cycles.
    max_cycles
        Positive intger, specifying the maximum number of micro cycles in a macro cycle.

    Returns
    -------
    dia_result
        Result of diabatization containing various quantities, e.g., the adiabatic-
        diabatic-transformation matrix.
    """
    assert nruns > 0
    assert df_tensor.ndim == 3
    nstates = df_tensor.shape[1]

    # Initial rotation matrices
    U0s = np.zeros((nruns, nstates, nstates))
    # Final, hopefully diabatized rotation matrices
    Us = np.zeros((nruns, nstates, nstates))
    jac_results = list()

    # Assume orthgonal states when no overlap matrix is given
    if overlap_matrix is None:
        overlap_matrix = np.eye(nstates)
    cost_funcs = np.zeros(nruns)
    diis_cycles = min(5, nstates)
    for i in range(nruns):
        # Generate random rotation matrix
        U0s[i] = get_random_U(nstates)
        # Use a copy of the rotation matrix as it is updated inplace
        jac_result = edmiston_ruedenberg(
            U0s[i].copy(),
            df_tensor,
            overlap_matrix,
            diis_cycles=diis_cycles,
            max_cycles=max_cycles,
        )
        jac_results.append(jac_result)
        if jac_result.is_converged:
            cost_funcs[i] = jac_result.P
            Us[i] = jac_result.C

    # Report result w/ highest value of cost_function
    max_ind = cost_funcs.argmax()
    best_jac_res = jac_results[max_ind]
    if not best_jac_res.is_converged:
        warnings.warn(
            f"Run {max_ind} has the maximum cost function value, but the optimization "
            f"did not converge after {max_cycles} cycles!"
        )
    return dia_result_from_jac_result(
        adia_ens,
        best_jac_res,
        L_tensor=df_tensor,
    )
