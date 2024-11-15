#!/usr/bin/env python

# [2] https://doi.org/10.1016/j.molliq.2021.115880
#     Individual tuning of solvent parameters -
#     from organic solvents to ionic liquids
#     Wankmüller, Berghold, Landgraf, 2021

import argparse
import sys

import jax
import jax.scipy as jsp
import jax.numpy as jnp
import numpy as np
import scipy as sp

from pysisyphus.constants import AU2EV, KBAU
from pysisyphus.helpers_pure import highlight_text
from pysisyphus.diabatization import logger
from pysisyphus.diabatization.coulomb import ERTemplate
from pysisyphus.diabatization.helpers import fmt_tensor
from pysisyphus.diabatization.results import (
    dia_result_from_jac_result,
    DiabatizationResult,
)
from pysisyphus.wavefunction.localization import JacobiSweepResult


__jac_func = None
__hess_func = None


def jac_func(*args, **kwargs):
    global __jac_func

    if __jac_func is None:
        __jac_func = jax.grad(er_eta)
    return __jac_func(*args, **kwargs)


def hess_func(*args, **kwargs):
    global __hess_func
    if __hess_func is None:
        __hess_func = jax.hessian(er_eta)
    return __hess_func(*args, **kwargs)


def rot_mat_from_skew_sym(k, n):
    """Skew-symmetric matrix of shape (n, n) from upper triangular values."""
    delta = jnp.zeros((n, n))
    triu = jnp.triu_indices(n, k=1)
    # Set supplied values at upper triangular part of the matrix.
    # The diagonal stays zero.
    delta = delta.at[triu].set(k)
    delta = delta - delta.T
    # Create rotation matrix from the exponential of a skew-symmetric matrix.
    U = jsp.linalg.expm(delta)
    return U


def er_eta(k, R, adia_ens, C, T=298.15):
    nstates = R.shape[0]
    U = rot_mat_from_skew_sym(k, nstates)
    # This contraction will be cheap, as nstates will be a small number.
    R_rot = jnp.einsum("Jj,Kk,Ll,Mm,JKLM->jklm", U, U, U, U, R, optimize="greedy")

    # Here it is important to use the energy differences between the supplied
    # electronic states and not the excitation energies from the ground state.
    adia_ens = adia_ens - adia_ens.min()
    beta = 1 / (KBAU * T)

    RAAAA = jnp.einsum("AAAA->A", R_rot)
    adia_ens_au = adia_ens / AU2EV
    EA = jnp.einsum("IA,I->A", U**2, adia_ens_au, optimize="greedy")
    exp_arg = -beta * (EA - C / 2 * RAAAA)
    cost_func = jnp.exp(exp_arg).sum()
    return cost_func


def rot_mat_from_er_eta(R, adia_ens, C, T):
    # TODO: revert this afterwards?!
    jax.config.update("jax_enable_x64", True)

    nstates = len(adia_ens)
    n = R.shape[0]
    nk = n * (n + 1) // 2 - n
    k0 = np.zeros(nk)
    k0 = np.random.rand(nk)

    eigval_thresh = 1e-8
    displ = 0.1
    k_cur = k0

    def fun(k):
        """As scipy only minimizes we return the negative of the function."""
        return -er_eta(k, R=R, adia_ens=adia_ens, C=C, T=T)

    max_cycles = 5
    logger.info("Edmiston-Ruedenberg-Eta cost-function maximization")
    for i in range(max_cycles):
        logger.info(f"Cost-function before macro cycle {i}: {-1*fun(k_cur): >12.6f}")
        res = sp.optimize.minimize(
            fun,
            x0=k_cur,
            method="Newton-CG",
            jac=lambda k: -jac_func(k, R=R, adia_ens=adia_ens, C=C, T=T),
            hess=lambda k: -hess_func(k, R=R, adia_ens=adia_ens, C=C, T=T),
            # callback=callback,
        )
        # logger.debug(res)
        k_cur = res.x

        hessian = hess_func(k_cur, R=R, adia_ens=adia_ens, C=C, T=T)
        w, v = np.linalg.eigh(hessian)
        non_neg = w > eigval_thresh
        logger.info(
            f"Macro cycle {i} converged in micro cycle {res.nit:>3d}. "
            f"There are {non_neg.sum()} positive eigenvalues."
        )
        is_converged = res.success and non_neg.sum() == 0
        if is_converged:
            logger.info("Macro cycles converged!")
            break
        # Displace current k along eigenvectors belonging to the positive eigenvalues
        non_neg_eigvecs = v[:, non_neg]
        for eigvec in non_neg_eigvecs.T:
            k_cur += displ * eigvec
    else:
        raise Exception(f"Diabatization did not converge in {max_cycles} cycles.")
    logger.info(f"Final cost-function: {-1*fun(k_cur): >12.6f}")

    # Final rotation matrix
    U = rot_mat_from_skew_sym(k_cur, n=nstates)

    msg = ERTemplate.render(
        name=highlight_text("Edmiston-Ruedenberg-Eta-diabatization"),
        R=R,
        fmt_tensor=fmt_tensor,
    )
    logger.info(msg)

    jac_res = JacobiSweepResult(
        is_converged=is_converged,
        cur_cycle=i,
        C=U,
        # Use actual cost function, not its negative
        P=-fun(k_cur),
    )
    return jac_res


def edmiston_ruedenberg_eta_diabatization(
    adia_ens: np.ndarray,
    R_tensor: np.ndarray,
    pekar: float,
    temperature: float,
) -> DiabatizationResult:
    # TODO: handle creation of these functions differently
    jac_res = rot_mat_from_er_eta(R_tensor, adia_ens, pekar, temperature)
    return dia_result_from_jac_result(
        kind="ereta",
        adia_ens=adia_ens,
        jac_res=jac_res,
    )


def callback(intermediate_result):
    # Report actual value
    act_val = -1 * intermediate_result.fun
    logger.info(f"{act_val: >12.6f}")


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("fn", type=str)
    parser.add_argument("T", default=298.15, type=float, help="Temperature in K.")
    parser.add_argument(
        "C",
        type=float,
        help=(
            "Pekar factor. Typical values range from 0.0 for non-polar "
            "solvents up to 0.55 for water at 25 °C."
        ),
    )
    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    fn = args.fn
    T = args.T
    C = args.C

    # Load data from ER-results
    data = np.load(fn)
    # Coulomb tensor
    R = data["R_tensor"]
    # Adiabatic excitation energies
    adia_ens = data["adia_ens"]

    U = rot_mat_from_er_eta(R, adia_ens, C, T)

    logger.info("Final rotation matrix U")
    with np.printoptions(precision=8, linewidth=240):
        logger.info(U)


if __name__ == "__main__":
    run()
