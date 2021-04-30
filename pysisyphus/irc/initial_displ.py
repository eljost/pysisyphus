# [1] https://pubs.acs.org/doi/10.1021/j100338a027
#     Koeski, Gordon, 1989
# [2] https://aip.scitation.org/doi/abs/10.1063/1.459634
#     Page, Doubleday, McIver, 1990

import numpy as np
from scipy.optimize import bisect


def get_curv_vec(H, Gv, v0, w0):
    v0 = v0[:, None]
    I = np.eye(v0.size)
    first = np.linalg.pinv(2 * w0 * I - H, rcond=1e-8)
    second = Gv.dot(v0) - v0.T.dot(Gv).dot(v0) * v0
    v1 = first.dot(second)
    return v1.flatten()


def taylor_closure(H, Gv, v0, v1, w0):
    """Taylor expansion of energy to 3rd order.

    dx(ds) = ds*v0 + ds**2*v1 / 2
    dE(ds) = dx^T H dx / 2 + dx^T [Gv] dx / 6

    H = Hessian
    Gv = 3rd derivative of energy along v0
    v0 = Transition vector
    v1 = Curvature vector
    w0 = Eigenvalue belonging to v0
    """

    # Precontract some values that will be reused
    v0v1 = v0.dot(v1)
    v1Hv1 = v1.dot(H).dot(v1)
    v0Gvv0 = v0.dot(Gv).dot(v0)
    v0Gvv1 = v0.dot(Gv).dot(v1)
    v1Gvv1 = v1.dot(Gv).dot(v1)

    def dE(ds):
        ds2 = ds ** 2
        ds3 = ds2 * ds
        ds4 = ds2 * ds2

        # Δx^T H Δx / 2
        quad = (w0 * ds2 + ds3 * w0 ** 2 * v0v1 + (ds4 * v1Hv1) / 4) / 2
        # Δx^T [Gv] Δx / 6
        cubic = (ds2 * v0Gvv0 + ds3 * v0Gvv1 + ds4 * v1Gvv1 / 4) / 6
        return quad + cubic

    return dE


def third_deriv_fd(geom, vec, ds=0.001):
    """Third derivative of the energy in direction 'vec'."""

    def get_H(geom, coords):
        H = geom.get_energy_and_cart_hessian_at(coords)["hessian"]
        H = geom.mass_weigh_hessian(H)
        # Only project for multi-atomic geometries.
        if geom.coords.size > 3:
            H = geom.eckart_projection(H)
        return H

    delta = ds * vec
    plus = geom.coords + delta
    minus = geom.coords - delta
    H_plus = get_H(geom, plus)
    H_minus = get_H(geom, minus)
    G_vec = (H_plus - H_minus) / (2 * ds)
    return G_vec


def cubic_displ(H, v0, w0, Gv, dE):
    v1 = get_curv_vec(H, Gv, v0, w0)
    E_taylor = taylor_closure(H, Gv, v0, v1, w0)

    def func(ds):
        return E_taylor(ds) - dE

    ds, rr = bisect(func, 0, 1, full_output=True)
    step = ds * v0 + ds ** 2 * v1 / 2
    return step


def cubic_displ_for_geom(geom, dE=-5e-4):
    H = geom.mw_hessian
    # Only project for multi-atomic geometries.
    if geom.coords.size > 3:
        H = geom.eckart_projection(H)
    w, v = np.linalg.eigh(H)
    w0 = w[0]
    v0 = v[:, 0]
    Gv = third_deriv_fd(geom, v0)
    return cubic_displ(H, v0, w0, Gv, dE=dE)
