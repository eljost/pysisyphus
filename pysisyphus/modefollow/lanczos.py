# [1] https://doi.org/10.1063/1.1809574
#     Comparison of methods for ﬁnding saddle points without
#     knowledge of the ﬁnal states
#     Olsen, Kroes, Henkelman, Arnaldsson, Jonsson, 2004

from collections.abc import Callable
import logging
from typing import Optional


import numpy as np
import scipy as sp


def lanczos(
    ncoords: int,
    grad_getter: Callable,
    dx: float = 5e-3,
    dl: float = 1e-2,
    guess: Optional[np.ndarray] = None,
    max_cycles: int = 25,
    reortho: bool = True,
    logger: Optional[logging.Logger] = None,
    fix_sign: bool = True,
) -> tuple[float, np.ndarray, bool]:
    """Lanczos method to determine smallest eigenvalue & -vector.

    See [1] for description of algorithm.

    Parameters
    ----------
    ncoords
        Dimensionality of the problem, e.g., number of rows/columns of the Hessian.
    grad_getter
        Function that calculates the gradient for a given displacement. It is NOT
        called with full coordinates, BUT with a displacement. Calculation of the
        full coordinates must be handled inside the function. See 'grad_getter()'
        in 'geom_lanczos()' below.
    dx
        Step size for finite differences calculation.
    dl
        Eigenvalue convergence threshold. See eq. (8) in [1].
        When abs((w_min - w_min_prev) / w_min_prev) < dl, convergence is signalled.
    guess
        Initial guess vector for the lowest eigenvector.
    max_cycles
        Maximum number of Lanczos cycles. Every cycle requires one gradient
        calculation.
    reortho
        Whether reorthogonalization is carried out.
    logger
        Logger object.
    fix_sign
        If true the first item of every eigenvector will be >= 0.0, e.g.

    Returns
    -------
    w_min
        Smallest eigenvector.
    eigenvector
        Normalized eigenvector belonging to the smallest eigenvalue.
    converged
        Whether the Lanczos iterations converged.
    """
    assert max_cycles > 0
    assert ncoords > 0

    def log(msg):
        if logger is not None:
            logger.debug(msg)

    log("Lanczos Algorithm")
    if guess is None:
        guess = np.random.rand(ncoords)
    else:
        guess = np.array(guess)
    assert guess.size == ncoords
    guess = guess / np.linalg.norm(guess)

    r_prev = guess
    beta_prev = np.linalg.norm(r_prev)
    q_prev = np.zeros_like(r_prev)

    alphas = np.empty(max_cycles)
    betas = np.empty(max_cycles)
    Q = np.zeros((ncoords, max_cycles))
    # Gradient at current coordinates; does not change throughout the iterations.
    grad_l = grad_getter(np.zeros(ncoords))
    w_min_prev = float("nan")
    for i in range(max_cycles):
        # Normalize q
        q = r_prev / beta_prev
        Q[:, i] = q
        # Approximate action of Hessian on q (u = Hq) by finite differences.
        # Eq. (9) and (10) in [1]
        # Gradient at perturbed coordinates
        grad_k = grad_getter(dx * q)
        u = (grad_k - grad_l) / dx
        # Residue
        r = u - beta_prev * q_prev
        alpha = q.dot(r)
        alphas[i] = alpha
        r -= alpha * q
        # Reorthogonalization of r against the present Lanczos vectors in Q
        if reortho:
            r -= Q.dot(Q.T.dot(r))
        beta = np.linalg.norm(r)
        betas[i] = beta

        """
        # Construct tri-diagonal matrix T
        size = i + 1
        T = np.zeros((size, size))
        diag_inds = np.diag_indices(size)
        T[diag_inds] = alphas[:size]
        if size > 1:
            for j, b in enumerate(betas[:i], 1):
                k = j - 1
                T[j, k] = b
                T[k, j] = b
        # Diagonalize T
        w, v = np.linalg.eigh(T)
        """

        # Earlier versions (up to 63ecd5d810a0f27563ce51a700625bc751220b5c) used
        # np.linalg.eigh and the "full" T-matrix. But as T is tridiagonal we can
        # use a specialized function for diagonalization.
        w, v = sp.linalg.eigh_tridiagonal(d=alphas[: i + 1], e=betas[:i])
        w_min = w[0]
        v_min = v[:, 0]

        # Form eigenvector belonging to the lowest eigenvalue from linear
        # combination of Lanczos vectors.
        eigenvector = (v_min[:, None] * Q.T[: i + 1]).sum(axis=0)
        eigenvector /= np.linalg.norm(eigenvector)
        dot = eigenvector.dot(guess)

        # Check eigenvalue convergence. Eq. (8) in [1]
        w_diff = w_min - w_min_prev
        w_quot = w_diff / w_min_prev
        converged = (i > 0) and abs(w_quot) < dl
        # Report current minimum eigenvalue, the quotient and the overlap between
        # current eigenvector and original guess vector.
        log(
            f"Cycle {i: >3d}: w_min={w_min: .6f}, {w_quot=: 12.4e}, v_min·guess={dot: >10.6f}"
        )
        if converged:
            log("Converged")
            break

        # Values for next cycle
        beta_prev = beta
        q_prev = q
        r_prev = r
        w_min_prev = w_min

    # Adapt sign, so the first element is >= 0 (positive).
    if fix_sign:
        negative = eigenvector[0] <= 0.0
        if negative:
            eigenvector = -eigenvector
    return w_min, eigenvector, converged


def geom_lanczos(geom, *args, **kwargs):
    """Wraps Lanczos algorithm for use with Geometry objects."""
    coords = geom.coords.copy()

    def grad_getter(displacement):
        results = geom.get_energy_and_forces_at(coords + displacement)
        return -results["forces"]

    ncoords = geom.coords.size
    return lanczos(ncoords, grad_getter, *args, **kwargs)
