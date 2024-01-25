# [1] https://doi.org/10.1063/1.469408
#     Efficient molecular numerical integration schemes
#     Treutler, Ahlrichs, 1995
# [2] https://doi.org/10.1080/00268979300100651
#     Quadrature schemes for integrals of density functional theory
#     Murray, Handy, Laming, 1992


import math
from typing import Callable

import numpy as np
import scipy as sp


LOG2 = math.log(2.0)


def treutler_m4_map(x: np.ndarray, a=1.0, alpha=0.6, atomic_radius=1.0):
    """Map finite inteval -1 <= x <= 1 to 0 <= r <= oo.

    For -1 <= x <= 1 a must be 1.0; for 0 <= x <= 1 a must be 0.0.

    Corresponds to scheme M4 in [1], eq. (19).

    Parameters
    ----------
    x
        Array of values in the interval [-1.0, 1.0].
    a
        Parameter a; see eq. (20) in [1]. Depends on the integration limits.
    alpha
        Parameter alpha with optimal value of 0.6 (see eq. (21) in [1]) in
        combination with Chebyshev 2nd kind roots and weights. As alpha is
        increased grid points are shifted towards r = 0 and r = oo.
    atomic_radius
        Additional parameter to better tailor the distribution of radial
        points.

    Returns
    -------
    r
        1d-array of r-values mapped from [-1, 1] to [0, oo].
    drdx
        1d-array of derivatives dr/dx that are required for the change of
        variables in the integration.
    """
    apx = a + x
    apx_alpha = apx**alpha
    log_ = np.log((a + 1.0) / (1.0 - x))
    r = atomic_radius * apx_alpha / LOG2 * log_
    xm1 = x - 1.0
    drdx = (
        atomic_radius
        * (alpha * apx_alpha * xm1 * log_ - apx_alpha * apx)
        / (apx * xm1 * LOG2)
    )
    return r, drdx


def chebyshev_2nd_kind(n: int, atomic_radius: float = 1.0):
    """Roots and weights using Chebyshev quadrature of the second kind.

    As outlined in [1].

    Parameters
    ----------
    n
        Number of radial grid points.
    atomic_radius
        Radius of the atom, for which the radial grid is generated.

    Returns
    -------
    roots_weight
        2d array of shape (n, 2) with the first column containing
        the radii and the second column containing the integration weights.
    """
    x, w = sp.special.roots_chebyu(n=n)
    inv_weight_func = 1 / np.sqrt(1 - x**2)
    r, drdx = treutler_m4_map(x, atomic_radius=atomic_radius)
    w *= drdx * inv_weight_func
    return np.stack((r, w), axis=1)


def euler_maclaurin_dma(n: int, m_r=2, atomic_radius: float = 1.0):
    """Roots and weights as utilized in Stone's GDMA.

    Adapated from the 'radial_grid(alpha)' subroutine found in'atomic_grids.f90'
    of Stone's GDMA code. Acutally a grid with one point less is returned to
    avoid division by 0.0.

    NOTE: In contrast to chebyshev_2nd_kind() above, the weights in the equation
    below already seem to include a factor of r**2. We remove it, so the function
    behaves similar to the other radial integration approaches.

    Parameters
    ----------
    n
        Positive integer > 1; number of points + 1.
    m_r
        Exponent, m_r = 2 was found to be sufficient. See Section 5.1 of [2]
        for a discussion.
    atomic_radius
        Positive float; corresponds to alpha in Stone's GDMA.

    Returns
    -------
    roots_weight
        2d array of shape (n, 2) with the first column containing
        the radii and the second column containing the integration weights.
    """
    tmp = np.arange(n - 1) + 1
    r = atomic_radius * (tmp / (n - tmp)) ** m_r
    prefact = m_r * n * atomic_radius**3
    w = prefact * tmp ** (3 * m_r - 1) / ((n - tmp) ** (3 * m_r + 1))
    # Remove r² factor to make this function consistent w/ the other approaches.
    w /= r**2
    return np.stack((r, w), axis=1)


def radint(
    func: Callable[[float], float],
    n: int,
    grid_func: Callable[[int], np.ndarray] = chebyshev_2nd_kind,
) -> float:
    """1d quadrature. Defaults to Chebyshev of the 2nd kind.

    Parameters
    ----------
    func
        Callable taking one argument.
    n
        Integer number of quadrature points.

    Returns
    -------
    res
        Floating point number holding the integration result.
    """

    rw = grid_func(n=n)
    r, w = rw.T
    return (func(r) * w).sum()
