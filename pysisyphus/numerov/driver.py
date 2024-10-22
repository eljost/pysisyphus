# [1] https://doi.org/10.1039/C6CP06698D
#     Pushing the limit for the grid-based treatment of Schrödinger's equation:
#     a sparse Numerov approach for one, two and three dimensional quantum problems
#     Kuenzer, Soraru, Hofer, 2016
# [2] https://doi.org/10.1119/1.4748813
#     Matrix Numerov method for solving Schrödinger’s equation
#     Pillai, Goglio, Walker, 2012
# [3] https://doi.org/10.1016/j.cplett.2019.04.016
#     A periodic Numerov approach applied to the torsional tunneling splitting
#     in hydrogen peroxide, aliphatic alcohols and phenol
#     Kuenzer, Hofer, 2019


from typing import Callable, Literal

import numpy as np
import scipy as sp


STENCILS = {
    2: np.array((1.0, -2, 1.0)),
    4: 1.0 / 12.0 * np.array((-1.0, 16.0, -30.0, 16.0, -1.0)),
    6: 1.0 / 180.0 * np.array((2.0, -27.0, 270.0, -490.0, 270.0, -27.0, 2.0)),
    8: (
        1.0
        / 5040.0
        * np.array(
            (-9.0, 128.0, -1008.0, 8064.0, -14350.0, 8064.0, -1008.0, 128.0, -9.0)
        )
    ),
    # There is a typo in the paper ... the central factor w/ subscript 0 must be 5269, not 5296!
    10: (
        1.0
        / 25200.0
        * np.array(
            [
                8.0,
                -125.0,
                1000.0,
                -6000.0,
                42000.0,
                -73766.0,
                42000.0,
                -6000.0,
                1000.0,
                -125.0,
                8.0,
            ]
        )
    ),
}


def off_diag_indices(n: int, k: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Return (off-)diagonal indices for a n-by-n matrix.

    To get the indices for a matrix of shape (500, 500) and a stencil size of 5
    this function must be called as 'off_diag_indices(500, 5)'.

    Parameters
    ----------
    n
        Length of the diagonal for k = 0 (number of rows/columns of the matrix).
    k
        Diagonal offset. For k = 0 indices for the main diagonal are returned.
        Positive k values 0 < k <= n - 1 yield indexs for diagonals above the main
        diagonal. Negative k value n - 1 <= k < 0 yield indexes for diagonals below
        the main diagonal.

    Returns
    -------
    diag_indices
        Tuple containing 2 1d-integer arrays indexing the desired diagonal.
    """

    abs_k = abs(k)
    ninds = n - abs_k
    assert ninds > 0, "abs(k) <= (n-1) must be fulfilled!"
    if k > 0:
        start = [0, abs_k]
    elif k < 0:
        start = [abs_k, 0]
    else:
        start = [0, 0]
    indices = np.zeros((2, ninds), dtype=int)
    ind = start
    for i in range(ninds):
        indices[:, i] = ind
        ind[0] += 1
        ind[1] += 1
    return indices[0], indices[1]


def get_A(n: int, d: float, coeffs: np.ndarray) -> np.ndarray:
    """Dense A-build.

    TODO: given the fact how simple the code for the sparse build
    is maybe we should just use it and convert the sparse matrix
    to a dense matrix, if desired. Then we could also delete
    off_diag_indices.

    Parameter
    ---------
    n
        Number of rows & columns in A.
    d
        Step size.
    stencil
        Finite-differences stencil.

    Returns
    -------
    A
        Dense A matrix of improved Numerov method.
    """
    coeffs = coeffs / d**2
    A = np.eye(n)
    start = (len(coeffs) - 1) // 2
    for k, od in enumerate(coeffs, -start):
        inds = off_diag_indices(n, k=k)
        A[inds] = od
    return A


def get_A_sparse(n: int, d: float, coeffs: np.ndarray) -> sp.sparse.dia_array:
    """Sparse A-build w/ DIAgonal storage.

    Parameter
    ---------
    n
        Number of rows & columns in A.
    d
        Step size.
    stencil
        Finite-differences stencil.
    periodic
        Whether periodic boundary calculations are used or not.

    Returns
    -------
    A
        Sparse A matrix of improved Numerov method with DIAgonal storage.
    """
    coeffs = coeffs / d**2
    start = (len(coeffs) - 1) // 2
    offsets = np.arange(-start, start + 1, dtype=int)
    data = np.repeat(coeffs[:, None], n, axis=1)
    return sp.sparse.dia_array((data, offsets), shape=(n, n))


def get_A_periodic(n: int, d: float, coeffs: np.ndarray) -> np.ndarray:
    """Dense A-build for a periodic potential.

    Parameter
    ---------
    n
        Number of rows & columns in A.
    d
        Step size.
    stencil
        Finite-differences stencil.

    Returns
    -------
    A
        Dense A matrix of improved Numerov method and a periodic potential.
    """
    ncoeffs = len(coeffs)
    assert ncoeffs % 2 == 1
    coeffs = coeffs / d**2
    mid = ncoeffs // 2
    stencil = np.arange(-mid, mid + 1, dtype=int)

    A = np.zeros((n, n))
    for i in range(n):
        stencil_i = stencil + i
        stencil_i[stencil_i >= n] %= n
        A[i, stencil_i] = coeffs
    return A


def get_A_periodic_sparse(n: int, d: float, coeffs: np.ndarray) -> sp.sparse.dia_array:
    """Sparse A-build for a periodic potential w/ DIAgonal storage.

    Parameter
    ---------
    n
        Number of rows & columns in A.
    d
        Step size.
    stencil
        Finite-differences stencil.
    periodic
        Whether periodic boundary calculations are used or not.

    Returns
    -------
    A
        Sparse A matrix of improved Numerov method with DIAgonal storage and a
        periodic potential.
    """
    A = sp.sparse.lil_array(get_A_periodic(n, d, coeffs))
    return A


def run(
    grid: np.ndarray,
    energy_getter: Callable[[int, float], float],
    mass: float,
    nstates: int = 50,
    accuracy: Literal[2, 4, 6, 8, 10] = 10,
    periodic: bool = False,
    normalize=True,
):
    """Improved 1d-Numerov method using a sparse matrix diagonalization.

    TODO: Normalize eigenvectors already here?

    Parameters
    ----------
    grid
        1d array containing an evenly spaced grid.
    energy_getter
        Callable that takes two arguments: an integer grid index and the floating point
        coordinate at this index. The index is useful in cases when all energies are
        already present and only have to be looked up.
    mass
        Mass in atomic units for the kinetic energy calculation.
    nstates
        Number of desired eigenstates to be calculated. Defaults to 50.
        As we use a sparse matrix diagonalization algorithm from scipy that converged
        the first nstats smallest eigenvalues nstates must be smaller than the number
        of grid points.
    accuracy
        Kind of stencil used for the improved Numerov-method. Must be one of
    periodic
        Boolean flag that indicates whether the potential is periodic or not.
    normalize
        Whether the returned eigenvectors are normalized so that |ψ|² integrates
        to 1.0. If set to False, than the Euclidian norm of the eigenvetors will
        be 1.0, but the integral over the grid is most likely != 1.0

    Returns
    -------
    w
        'nstates' smallest eigenvalues for given potential.
    v
        'nstates' eigenvectors belonging to the calculated eigenvectors.
    """
    d = grid[1] - grid[0]
    assert d > 0.0
    assert mass > 0.0
    stencil = STENCILS[accuracy]

    # Number of grid points
    n = len(grid)

    get_A_func = get_A_periodic_sparse if periodic else get_A_sparse
    Asp = get_A_func(n, d, stencil)
    # Calculate potential energy on grid
    V_dia = np.zeros_like(grid)
    for i in range(n):
        V_dia[i] = energy_getter(i, grid[i])
    Vsp = sp.sparse.dia_array((V_dia[None, :], [0]), shape=(n, n))
    # Build sparse Hamiltonian
    Hsp = -Asp / (2.0 * mass) + Vsp
    # Diagonalize sparse Hamiltonian to calculate the first 'nstates' smallest
    # eigenvalues and eigenvectors.
    w, v = sp.sparse.linalg.eigsh(Hsp, k=nstates, which="SM")

    if normalize:
        # In principle it should be enough to determine the normalization
        # factor for the first state and to apply it to all others, as it should be
        # entirely grid-size dependent.
        #
        # But given the fact that nstates will usually be a small number (nstates << 100)
        # it should be cheap to numerically intgrate and normalize each state on its own.
        for i, psi in enumerate(v.T):
            abs_psi2 = np.abs(psi) ** 2
            norm = sp.integrate.simpson(abs_psi2, x=grid, dx=d)
            psi /= np.sqrt(norm)
            v[:, i] = psi
    return w, v
