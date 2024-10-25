from collections.abc import Sequence
import functools
from typing import Callable, Literal, Optional

import numpy as np

from pysisyphus.executors import CloudpickleProcessPoolExecutor


def finite_difference_gradient(
    coords: np.ndarray,
    scalar_func: Callable,
    step_size: float = 1e-2,
):
    """Stripped down version of finite_difference_hessian.

    Both functions could probably be unified."""
    size = coords.size
    fd_gradient = np.zeros(size)
    zero_step = np.zeros(size)

    coeffs = ((-0.5, -1), (0.5, 1))
    for i, _ in enumerate(coords):
        step = zero_step.copy()
        step[i] = step_size

        def get_scalar(factor, displ):
            displ_coords = coords + step * displ
            scalar = scalar_func(displ_coords)
            return factor * scalar

        plus, minus = [get_scalar(factor, displ) for factor, displ in coeffs]
        fd_gradient[i] = (plus + minus) / step_size

    return fd_gradient


def finite_difference_hessian(
    coords: np.ndarray,
    grad_func: Callable[[np.ndarray], np.ndarray],
    step_size: float = 1e-2,
    acc: Literal[2, 4] = 2,
    callback: Optional[Callable] = None,
) -> np.ndarray:
    """Numerical Hessian from central finite gradient differences.

    See central differences in
      https://en.wikipedia.org/wiki/Finite_difference_coefficient
    for the different accuracies.
    """
    if callback is None:

        def callback(*args):
            pass

    accuracies = {
        # key: ((factor0, displacement_factor0), (factor1, displacement_factor1), ...)
        2: ((-0.5, -1), (0.5, 1)),  # 2 calculations
        4: ((1 / 12, -2), (-2 / 3, -1), (2 / 3, 1), (-1 / 12, 2)),  # 4 calculations
    }
    accs_avail = list(accuracies.keys())
    assert acc in accs_avail

    size = coords.size
    fd_hessian = np.zeros((size, size))
    zero_step = np.zeros(size)

    coeffs = accuracies[acc]
    for i, _ in enumerate(coords):
        step = zero_step.copy()
        # Change one coordinate at a time
        step[i] = step_size

        def get_grad(factor, displ):
            displ_coords = coords + step * displ
            callback(i)
            grad = grad_func(displ_coords)
            return factor * grad

        # Depending on the chosen accuracy, a different number of calculations
        # have to be carried out.
        grads = [get_grad(factor, displ) for factor, displ in coeffs]
        # Even though we are interested in the second derivative of the energy w.r.t.
        # the coordinates, we calculate it from the first derivative of the gradient.
        # So we only divide by step_size**1, not by step_size**2.
        fd = np.sum(grads, axis=0) / step_size
        fd_hessian[i] = fd

    # Symmetrize
    fd_hessian = (fd_hessian + fd_hessian.T) / 2
    return fd_hessian


def grad_func(index, coords, factor, calc, atoms, prepare_kwargs: dict):
    if index is not None:
        # calc.base_name = f"{calc.base_name}_num_hess"
        calc.calc_counter = index
    results = calc.get_forces(atoms, coords, **prepare_kwargs)
    gradient = -results["forces"]
    return factor * gradient


def finite_difference_hessian_mp(
    atoms: tuple[str],
    coords: np.ndarray,
    calc,
    step_size: float = 1e-2,
    acc: Literal[2, 4] = 2,
    serial: bool = False,
    prepare_kwargs: Optional[dict] = None,
) -> np.ndarray:
    """Finite-differences Hessian w/ multiprocessing support.

    See central differences in
      https://en.wikipedia.org/wiki/Finite_difference_coefficient
    for the different accuracies.
    """
    if prepare_kwargs is None:
        prepare_kwargs = {}

    try:
        org_pal = calc.pal
    except AttributeError:
        org_pal = 1
    # When multiprocessing is used all calculations will be carried out
    # on one core each. In serial calculations we'll restore the original
    # pal value.
    calc.pal = 1

    # Pre-supply some arguments so we later only have to supply the arguments
    # speecific to the displacments.
    gf = functools.partial(
        grad_func, calc=calc, atoms=atoms, prepare_kwargs=prepare_kwargs
    )

    accuracies = {
        # key: ((factor0, displacement_factor0), (factor1, displacement_factor1), ...)
        2: ((-0.5, -1), (0.5, 1)),  # 2 calculations
        4: ((1 / 12, -2), (-2 / 3, -1), (2 / 3, 1), (-1 / 12, 2)),  # 4 calculations
    }
    accs_avail = list(accuracies.keys())
    assert acc in accs_avail
    coeffs = accuracies[acc]

    size = coords.size
    fd_hessian = np.zeros((size, size))

    zero_step = np.zeros(size)
    row_inds = list()
    factors = list()
    all_displ_coords = list()
    # Prepare lists of displacments and associated factors, we can use them later
    # in map-like function calls.
    for i, _ in enumerate(coords):
        for factor, displ in coeffs:
            row_inds.append(i)
            factors.append(factor)
            # Change one coordinate at a time
            step = zero_step.copy()
            step[i] = displ * step_size
            all_displ_coords.append(coords + step)
    ndispls = len(all_displ_coords)

    if serial or (org_pal == 1):
        # In serial mode we dont have to modify calc_counter, as we operate on the
        # same object throughout.
        calc.pal = org_pal
        ndispls = [None] * ndispls
        grads = [
            gf(i, displ_coords, factor)
            for i, displ_coords, factor in zip(ndispls, all_displ_coords, factors)
        ]
    else:
        with CloudpickleProcessPoolExecutor(max_workers=org_pal) as executor:
            grads = list(executor.map(gf, range(ndispls), all_displ_coords, factors))
        # Update calc_counter, as the actual calculations were carried out on copies
        # of the calculator.
        calc.calc_counter += ndispls
    # Restore original pal value
    calc.pal = org_pal

    # Distribute the gradients on the appropriate Hessian rows
    for row_ind, grad in zip(row_inds, grads):
        fd_hessian[row_ind] += grad

    # Alternatively, the factors could already be divided by the step_size.
    fd_hessian /= step_size
    # Symmetrize
    fd_hessian = (fd_hessian + fd_hessian.T) / 2
    return fd_hessian


def periodic_finite_difference(
    x_ind: int, y: np.ndarray, dx: float, degree: int, coeffs: Sequence[float]
):
    """Finite differences for data on a periodic grid.

    Parameters
    ----------
    x_ind
        Non-negative integer indexing the point on the grid, where the derivative
        is evaluated.
    y
        Data array on a periodic evenely spaced grid.
    dx
        Grid spacing.
    degree
        Degree of the derivative. Positive integer >= 1.
    coeffs
        Sequence of floating point numbers containing Floating pi
    Returns
    -------
    """
    assert x_ind >= 0
    assert degree >= 1, "Degree must be >= 1 but got {degree}!"
    assert len(coeffs) % 2 == 1, "Length of coefficient sequence must be odd!"
    mid = (len(coeffs) - 1) // 2
    npoints = len(y)
    # The negative indices will wrap around.
    stencil = np.arange(-mid, mid + 1) + x_ind
    mask = stencil > npoints - 1
    stencil[mask] = stencil[mask] - npoints
    deriv = (y[stencil] * coeffs).sum() / dx**degree
    return deriv


def periodic_fd_2_8(x_ind: int, y: np.ndarray, dx: float):
    """2nd derivative for 8th-order accuracy from data on a periodic grid.

    For more information see the docstring of 'periodic_finite_difference'.
    """
    return periodic_finite_difference(
        x_ind,
        y,
        dx,
        degree=2,
        coeffs=[
            -1.0 / 560.0,
            8.0 / 315.0,
            -1.0 / 5.0,
            8.0 / 5.0,
            -205.0 / 72.0,
            8.0 / 5.0,
            -1.0 / 5.0,
            8.0 / 315.0,
            -1.0 / 560.0,
        ],
    )
