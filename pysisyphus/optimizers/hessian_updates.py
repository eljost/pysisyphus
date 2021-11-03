# [1] https://link.springer.com/article/10.1007/s00214-016-1847-3
#     Birkholz, 2016
# [2] Geometry optimization in Cartesian coordinates: Constrained optimization
#     Baker, 1992
# [3] https://epubs.siam.org/doi/pdf/10.1137/S1052623496306450
#     BFGS WITH UPDATE SKIPPING AND VARYING MEMORY
#     Kolda, 1998
# [4] https://link.springer.com/article/10.1186/1029-242X-2012-241
#     New cautious BFGS algorithm based on modified Armijo-type line search
#     Wan, 2012
# [5] Numerical optimization, 2nd ed.
#     Nocedal, Wright
# [6] https://arxiv.org/abs/2006.08877
#     Goldfarb, 2020
# [7] https://pubs.acs.org/doi/10.1021/acs.jctc.9b00869
#     Hermes, Zádor, 2019
# [8] https://doi.org/10.1002/(SICI)1096-987X(199802)19:3<349::AID-JCC8>3.0.CO;2-T
#     Bofill, 1998
# [9] http://dx.doi.org/10.1016/S0166-1280(02)00209-9
#     Bungay, Poirier


import numpy as np

from pysisyphus.optimizers.closures import bfgs_multiply


def bfgs_update(H, dx, dg):
    first_term = np.outer(dg, dg) / dg.dot(dx)
    second_term = H.dot(np.outer(dx, dx)).dot(H) / dx.dot(H).dot(dx)
    return first_term - second_term, "BFGS"


def damped_bfgs_update(H, dx, dg):
    """See [5]"""
    dxdg = dx.dot(dg)
    dxHdx = dx.dot(H).dot(dx)
    theta = 1
    if dxdg < 0.2 * dxHdx:
        theta = 0.8 * dxHdx / (dxHdx - dxdg)
    r = theta * dg + (1 - theta) * H.dot(dx)

    first_term = np.outer(r, r) / r.dot(dx)
    second_term = H.dot(np.outer(dx, dx)).dot(H) / dxHdx
    return first_term - second_term, "damped BFGS"


def double_damp(
    s, y, H=None, s_list=None, y_list=None, mu_1=0.2, mu_2=0.2, logger=None
):
    """Double damped step 's' and gradient differences 'y'.

    H is the inverse Hessian!
    See [6]. Potentially updates s and y. y is only
    updated if mu_2 is not None.

    Parameters
    ----------
    s : np.array, shape (N, ), floats
        Coordiante differences/step.
    y : np.array, shape (N, ), floats
        Gradient differences
    H : np.array, shape (N, N), floats, optional
        Inverse Hessian.
    s_list : list of nd.array, shape (K, N), optional
        List of K previous steps. If no H is supplied and prev_ys is given
        the matrix-vector product Hy will be calculated through the
        two-loop LBFGS-recursion.
    y_list : list of nd.array, shape (K, N), optional
        List of K previous gradient differences. See s_list.
    mu_1 : float, optional
        Parameter for 's' damping.
    mu_2 : float, optional
        Parameter for 'y' damping.
    logger : logging.Logger, optional
        Logger to be used.

    Returns
    -------
    s : np.array, shape (N, ), floats
        Damped coordiante differences/step.
    y : np.array, shape (N, ), floats
        Damped gradient differences
    """
    sy = s.dot(y)
    # Calculate Hy directly
    if H is not None:
        Hy = H.dot(y)
    # Calculate Hy via BFGS_multiply as in LBFGS
    else:
        Hy = bfgs_multiply(s_list, y_list, y, logger=logger)
    yHy = y.dot(Hy)

    theta_1 = 1
    damped_s = ""
    if sy < mu_1 * yHy:
        theta_1 = (1 - mu_1) * yHy / (yHy - sy)
        s = theta_1 * s + (1 - theta_1) * Hy
        if theta_1 < 1.0:
            damped_s = ", damped s"
    msg = f"damped BFGS\n\ttheta_1={theta_1:.4f} {damped_s}"

    # Double damping
    damped_y = ""
    if mu_2 is not None:
        sy = s.dot(y)
        ss = s.dot(s)
        theta_2 = 1
        if sy < mu_2 * ss:
            theta_2 = (1 - mu_2) * ss / (ss - sy)
        y = theta_2 * y + (1 - theta_2) * s
        if theta_2 < 1.0:
            damped_y = ", damped y"
        msg = "double " + msg + f"\n\ttheta_2={theta_2:.4f} {damped_y}"

    if logger is not None:
        logger.debug(msg.capitalize())

    return s, y


def sr1_update(z, dx):
    return np.outer(z, z) / z.dot(dx), "SR1"


def psb_update(z, dx):
    first_term = (np.outer(dx, z) + np.outer(z, dx)) / dx.dot(dx)
    sec_term = dx.dot(z) * np.outer(dx, dx) / dx.dot(dx) ** 2
    return first_term - sec_term, "PSB"


def flowchart_update(H, dx, dg):
    # See [1], Sec. 2, equations 1 to 3
    z = dg - H.dot(dx)
    sr1_quot = z.dot(dx) / (np.linalg.norm(z) * np.linalg.norm(dx))
    bfgs_quot = dg.dot(dx) / (np.linalg.norm(dg) * np.linalg.norm(dx))
    if sr1_quot < -0.1:
        update, key = sr1_update(z, dx)
    elif bfgs_quot > 0.1:
        update, key = bfgs_update(H, dx, dg)
    else:
        update, key = psb_update(z, dx)
    return update, key


def mod_flowchart_update(H, dx, dg):
    # This version seems to work too ... at least for minimizations
    # starting from a geometry near a transition state. Interesing.
    z = dg - H.dot(dx)
    quot = z.dot(dx) / (np.linalg.norm(z) * np.linalg.norm(dx))
    if quot < -0.1:
        update, key = bfgs_update(H, dx, dg)
    elif quot > 0.1:
        update, key = sr1_update(z, dx)
    else:
        update, key = psb_update(z, dx)
    return update, key


def bofill_update(H, dx, dg):
    z = dg - H.dot(dx)

    # Symmetric, rank-one (SR1) update
    sr1, _ = sr1_update(z, dx)

    # Powell (PSB) update
    powell, _ = psb_update(z, dx)

    # Bofill mixing-factor
    mix = z.dot(dx) ** 2 / (z.dot(z) * dx.dot(dx))

    # Bofill update
    bofill_update = (mix * sr1) + (1 - mix) * (powell)

    return bofill_update, "Bofill"


def ts_bfgs_update(H, dx, dg):
    """As described in [7]"""
    dx = dx[:, None]
    dg = dg[:, None]
    j = dg - H @ dx
    jdx = j.T @ dx
    # Diagonalize Hessian, to construct positive definite version of it
    w, v = np.linalg.eigh(H)
    Hdx = np.abs(w) * v @ (v.T @ dx)
    M = dg @ dg.T + Hdx @ Hdx.T
    dxTM = dx.T @ M
    u = np.linalg.solve(dxTM @ dx, dxTM).T
    juT = j @ u.T
    ts_bfgs_update = juT + juT.T - jdx * u @ u.T
    return ts_bfgs_update, "TS-BFGS"


def ts_bfgs_update_org(H, dx, dg):
    """Do not use! Implemented as described in the 1998 bofill paper [8].

    This does not seem to work too well."""
    dx = dx[:, None]
    dg = dg[:, None]
    u = dg
    j = H @ dx
    j = u - j
    # jTdx = float(j.T @ dx)
    jTdx = j.T @ dx
    # dxTdx = float(dx.T @ dx)
    dxTdx = dx.T @ dx
    # jTj = float(j.T @ j)
    jTj = j.T @ j
    phi = jTdx ** 2 / dxTdx / jTj
    u = phi * dg * dg.T @ dx
    w, v = np.linalg.eigh(H)
    u = u + (1 - phi) * np.abs(w) * v @ (v.T @ dx)
    u = u / (u.T @ dx)
    juT = j @ u.T
    ts_bfgs_update = juT + juT.T - jTdx * u @ u.T
    return ts_bfgs_update, "TS-BFGS"


def ts_bfgs_update_revised(H, dx, dg):
    """TS-BFGS update as described in [9].

    Better than the original formula of Bofill, worse than the implementation
    in [7]. a is caluclated as described in the footnote 1 on page 38. Eq. (8)
    looks suspicious as it contains the inverse of a vector?! As also outlined
    in the paper abs(a) is used (|a| in the paper)."""

    dx = dx[:, None]
    dg = dg[:, None]
    dgTdg = dg.T @ dg
    dgTdx = dg.T @ dx
    a = (dgTdg - dg.T @ H @ dx) / (dgTdg * dgTdx)
    a = abs(a)

    # Diagonalize Hessian, to construct positive definite version of it
    w, v = np.linalg.eigh(H)
    H_pos_dx = np.abs(w) * v @ (v.T @ dx)
    # Mixing factor
    j = dg - H @ dx
    jTdx = j.T @ dx
    dxTdx = dx.T @ dx
    jTj = j.T @ j
    phi = jTdx ** 2 / dxTdx / jTj
    u = ((1 - phi) * H_pos_dx + a * phi * dgTdx * dg) / (
        (1 - phi) * dx.T @ H_pos_dx + phi * a * dgTdx ** 2
    )

    juT = j @ u.T
    ts_bfgs_update = juT + juT.T - jTdx * u @ u.T
    return ts_bfgs_update, "TS-BFGS"
