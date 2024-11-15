import itertools as it
from typing import Optional
import warnings

from jinja2 import Template
import numpy as np

from pysisyphus.diabatization.helpers import fmt_tensor, get_random_U
from pysisyphus.diabatization.results import (
    DiabatizationResult,
    dia_result_from_jac_result,
)
from pysisyphus.helpers_pure import highlight_text
from pysisyphus.diabatization import logger
from pysisyphus.wavefunction.localization import JacobiSweepResult, edmiston_ruedenberg


ERTemplate = Template(
    """{{ name }}

Coulomb-Tensor R
----------------
{{ fmt_tensor(R) }}

"""
)


def edmiston_ruedenberg_jacobi_sweeps(
    R: np.ndarray,
    random: bool = False,
    U=None,
    max_cycles: Optional[int] = None,
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
        fmt_tensor=fmt_tensor,
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
        logger.info(f"Macro cycle {i: >3d}: {D=: >14.6f}, {dD=: >+12.6e}")
        if converged := (abs(dD) <= dP_thresh):
            logger.info("Converged!")
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
            # logger.info(
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
            logger.warning("A**2 + B**2 is very small. Indicating convergence.")
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
        # logger.info(f"\tRotation of {j=} and {k=} by {rad=: 10.6e} rad ({deg: >10.6f}°).")
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
    kind = "Edmiston-Ruedenberg"
    return dia_result_from_jac_result(
        kind,
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
    kind = "Edmiston-Ruedenberg"
    return dia_result_from_jac_result(
        kind,
        adia_ens,
        best_jac_res,
        L_tensor=df_tensor,
    )
