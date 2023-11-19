import distributed
import numpy as np
from typing import Callable, Literal, Optional


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


def finite_difference_hessian_pal(
    coords: np.ndarray,
    grad_func: Callable[[np.ndarray], np.ndarray],
    step_size: float = 1e-2,
    acc: Literal[2, 4] = 2,
    client: distributed.Client = None,
) -> np.ndarray:
    """Numerical Hessian from central finite gradient differences.

    See central differences in
      https://en.wikipedia.org/wiki/Finite_difference_coefficient
    for the different accuracies.
    """

    def get_grad(calc_number, step, factor):
        displ_coords = coords + step
        grad = grad_func(calc_number, displ_coords)
        return factor * grad

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
    steps = list()
    for i, _ in enumerate(coords):
        for factor, displ in coeffs:
            row_inds.append(i)
            factors.append(factor)
            # Change one coordinate at a time
            step = zero_step.copy()
            step[i] = displ * step_size
            steps.append(step)
    calc_numbers = list(range(len(steps)))

    if client is not None:
        futures = client.map(get_grad, calc_numbers, steps, factors, pure=False)
        grads = client.gather(futures)
    else:
        grads = [
            get_grad(calc_number, step, factor)
            for calc_number, step, factor in zip(calc_numbers, steps, factors)
        ]

    for row_ind, grad in zip(row_inds, grads):
        fd_hessian[row_ind] += grad

    # Alternatively, the factors could already be divided by the step_size.
    fd_hessian /= step_size
    # Symmetrize
    fd_hessian = (fd_hessian + fd_hessian.T) / 2
    return fd_hessian
