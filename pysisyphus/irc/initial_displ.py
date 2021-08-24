# [1] https://pubs.acs.org/doi/10.1021/j100338a027
#     Koeski, Gordon, 1989
# [2] https://aip.scitation.org/doi/abs/10.1063/1.459634
#     Page, Doubleday, McIver, 1990

from collections import namedtuple

import h5py
import numpy as np
from scipy.optimize import bisect

from pysisyphus.constants import AU2KJPERMOL


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


ThirdDerivResult = namedtuple(
    "ThirdDerivResult",
    "coords_plus energy_plus H_plus coords_minus energy_minus H_minus G_vec vec ds",
)


def third_deriv_fd(geom, vec, ds=0.001):
    """Third derivative of the energy in direction 'vec'."""

    def get_H(geom, coords):
        results = geom.get_energy_and_cart_hessian_at(coords)
        energy = results["energy"]
        H = results["hessian"]
        H_mw = geom.mass_weigh_hessian(H)
        # Only project for multi-atomic geometries.
        if geom.coords.size > 3:
            H_mw = geom.eckart_projection(H_mw)
        return H_mw, H, energy

    delta = ds * vec
    plus = geom.coords + delta
    minus = geom.coords - delta
    H_mw_plus, H_plus, energy_plus = get_H(geom, plus)
    H_mw_minus, H_minus, energy_minus = get_H(geom, minus)
    G_vec = (H_mw_plus - H_mw_minus) / (2 * ds)

    third_deriv_res = ThirdDerivResult(
        coords_plus=plus,
        energy_plus=energy_plus,
        H_plus=H_plus,
        coords_minus=minus,
        energy_minus=energy_minus,
        H_minus=H_minus,
        G_vec=G_vec,
        vec=vec,
        ds=ds,
    )
    return G_vec, third_deriv_res


def cubic_displ(H, v0, w0, Gv, dE):
    """
    According to Eq. (26) in [2] v1 does not depend on the sign of v0.
        v1 = (F0 - 2v0^T F0 v0 I)⁻¹ x ([G0v0] - v0^T [G0v0] v0 I) v0
    The first term is obviously independent of v0's sign. Using v0' = -v0 the
    second term becomes
        ([G0v0'] - v0'^T [G0v0'] v0' I) v0'
        (-[G0v0] - v0^T [G0v0'] v0 I) v0'
        (-[G0v0] + v0^T [G0v0] v0 I) v0'
        -(-[G0v0] + v0^T [G0v0] v0 I) v0
        ([G0v0] - v0^T [G0v0] v0 I) v0
    Strictly speaking the contraction [G0v0] should depend on the sign of v0.
    In the current implementation, the direction of v0 is not taken into account,
    but
        get_curv_vec(H, Gv, v0, w0) == get_curv_vec(H, -Gv, -v0, w0) .
    But somehow the Taylor expansion gives bogus results when called with -Gv and -v0...
    """

    assert dE < 0.0, f"Supplied dE={dE:.6f} is positive but it must be negative!"
    assert w0 < 0.0, f"Expected first eigenvalue to be negative but it is w0={w0:.6e}!"

    v1 = get_curv_vec(H, Gv, v0, w0)
    E_taylor = taylor_closure(H, Gv, v0, v1, w0)

    def func(ds):
        return E_taylor(ds) - dE

    def prepare_bisect(x0, theta=1.25, max_cycles=20):
        """Determine (lower) upper bound for scipy.optimize.bisect."""
        assert theta > 1.0

        ds = x0
        dE_min = func(0.0)
        dE_prev = dE_min
        ds_min = ds
        # Grow until we find an upper (lower) bound of the interval
        for _ in range(max_cycles):
            dE = func(ds)
            # print(
                # f"ds={ds:.4f} dE={dE*AU2KJPERMOL:.3f} dE_min={dE_min*AU2KJPERMOL:.3f}"
            # )

            if dE <= 0.0:
                break
            # Keep best guess
            elif dE <= dE_min:
                dE_min = dE
                ds_min = ds
            # Return (yet) best guess when function value grows again
            elif dE >= dE_prev:
                ds = ds_min
                break

            dE_prev = dE
            ds *= theta  # Grow ds
        return ds

    def bisect_(ds0):
        try:
            ds_, rr = bisect(func, a=0.0, b=ds0, full_output=True)
        # Will be raised when f(a) and f(b) have the same sign.
        except ValueError:
            ds_ = ds0
        return ds_

    plus_bound = prepare_bisect(0.1)
    ds_plus = bisect_(plus_bound)

    minus_bound = prepare_bisect(-0.1)
    ds_minus = bisect_(minus_bound)

    # import matplotlib.pyplot as plt
    # dss = np.linspace(-1, 1, num=51)
    # # dss = np.linspace(-6, 6, num=200)
    # E0 = E_taylor(0.0)
    # Es = E_taylor(dss) - E0
    # Es *= AU2KJPERMOL
    # Emp = (np.array((E_taylor(ds_minus), E_taylor(ds_plus))) - E0) * AU2KJPERMOL
    # _, ax = plt.subplots()
    # ax.plot(dss, Es, "o-")
    # ax.axvline(0.0, c="k", ls="--")
    # ax.axhline(dE*AU2KJPERMOL, c="k", ls=":")
    # ax.scatter((ds_minus, ds_plus), Emp, s=75, marker="x", c="r", zorder=3)
    # ax.set_xlabel("ds")
    # ax.set_ylabel("dE / kJ mol⁻¹")
    # plt.show()

    def step(ds):
        return ds * v0 + ds ** 2 * v1 / 2

    step_plus = step(ds_plus)
    step_minus = step(ds_minus)
    return step_plus, step_minus


def cubic_displ_for_h5(h5_fn="third_deriv.h5", dE=-5e-4):
    with h5py.File(h5_fn, "r") as handle:
        H_mw = handle["H_mw"][:]
        Gv = handle["G_vec"][:]

    w, v = np.linalg.eigh(H_mw)
    w0 = w[0]
    v0 = v[:, 0]
    return cubic_displ(H_mw, v0, w0, Gv, dE)


def cubic_displ_for_geom(geom, dE=-5e-4):
    H = geom.mw_hessian
    # Only project for multi-atomic geometries.
    if geom.coords.size > 3:
        H = geom.eckart_projection(H)
    w, v = np.linalg.eigh(H)
    # Transition vector (imaginary mode) and corresponding eigenvalue
    v0 = v[:, 0]
    w0 = w[0]
    Gv, third_deriv_res = third_deriv_fd(geom, v0)
    step_plus, step_minus = cubic_displ(H, v0, w0, Gv, dE=dE)
    return step_plus, step_minus, third_deriv_res
