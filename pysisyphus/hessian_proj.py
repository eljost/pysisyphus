# [1] https://doi.org/10.1063/1.468630
#     Reaction‐path potential and vibrational frequencies in terms
#     of curvilinear internal coordinates
#     Jackels, Gu, Truhlar, 1995
# [2] https://doi.org/10.1021/jp9724028
#     Reaction-Path Dynamics in Redundant Internal Coordinates
#     Chuang, Truhlar, 1998
# [3] https://onlinelibrary.wiley.com/doi/full/10.1002/jcc.10089
#     Quantum chemical calculation of vibrational spectra of large
#     molecules — Raman and IR spectra for Buckminsterfullerene
#     Neugebauer, Reiher, Hess, 2002


import warnings

import numpy as np

from pysisyphus.helpers_pure import rms


def inertia_tensor(coords3d, masses):
    """Inertita tensor.

                          | x² xy xz |
    (x y z)^T . (x y z) = | xy y² yz |
                          | xz yz z² |
    """
    x, y, z = coords3d.T
    squares = np.sum(coords3d**2 * masses[:, None], axis=0)
    I_xx = squares[1] + squares[2]
    I_yy = squares[0] + squares[2]
    I_zz = squares[0] + squares[1]
    I_xy = -np.sum(masses * x * y)
    I_xz = -np.sum(masses * x * z)
    I_yz = -np.sum(masses * y * z)
    I = np.array(((I_xx, I_xy, I_xz), (I_xy, I_yy, I_yz), (I_xz, I_yz, I_zz)))
    return I


def get_trans_rot_vectors(cart_coords, masses, rot_thresh=1e-6):
    """Vectors describing translation and rotation.

    These vectors are used for the Eckart projection by constructing
    a projector from them.

    See Martin J. Field - A Pratcial Introduction to the simulation
    of Molecular Systems, 2007, Cambridge University Press, Eq. (8.23),
    (8.24) and (8.26) for the actual projection.

    See also https://chemistry.stackexchange.com/a/74923.

    Parameters
    ----------
    cart_coords : np.array, 1d, shape (3 * atoms.size, )
        Atomic masses in amu.
    masses : iterable, 1d, shape (atoms.size, )
        Atomic masses in amu.

    Returns
    -------
    ortho_vecs : np.array(6, 3*atoms.size)
        2d array containing row vectors describing translations
        and rotations.
    """

    coords3d = np.reshape(cart_coords, (-1, 3))
    total_mass = masses.sum()
    com = 1 / total_mass * np.sum(coords3d * masses[:, None], axis=0)
    coords3d_centered = coords3d - com[None, :]

    I = inertia_tensor(coords3d, masses)
    _, Iv = np.linalg.eigh(I)
    Iv = Iv.T

    masses_rep = np.repeat(masses, 3)
    sqrt_masses = np.sqrt(masses_rep)
    num = len(masses)

    def get_trans_vecs():
        """Mass-weighted unit vectors of the three cartesian axes."""

        for vec in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            _ = sqrt_masses * np.tile(vec, num)
            yield _ / np.linalg.norm(_)

    def get_rot_vecs():
        """As done in geomeTRIC."""

        rot_vecs = np.zeros((3, cart_coords.size))
        # p_vecs = Iv.dot(coords3d_centered.T).T
        for i in range(masses.size):
            p_vec = Iv.dot(coords3d_centered[i])
            for ix in range(3):
                rot_vecs[0, 3 * i + ix] = Iv[2, ix] * p_vec[1] - Iv[1, ix] * p_vec[2]
                rot_vecs[1, 3 * i + ix] = Iv[2, ix] * p_vec[0] - Iv[0, ix] * p_vec[2]
                rot_vecs[2, 3 * i + ix] = Iv[0, ix] * p_vec[1] - Iv[1, ix] * p_vec[0]
        rot_vecs *= sqrt_masses[None, :]
        return rot_vecs

    trans_vecs = list(get_trans_vecs())
    rot_vecs = np.array(get_rot_vecs())
    # Drop vectors with vanishing norms
    rot_vecs = rot_vecs[np.linalg.norm(rot_vecs, axis=1) > rot_thresh]
    tr_vecs = np.concatenate((trans_vecs, rot_vecs), axis=0)
    tr_vecs = np.linalg.qr(tr_vecs.T)[0].T
    return tr_vecs


def get_projector(vecs):
    """Vectors must be in rows.

    # P = I - sum_i vec_i @ vec_i.T
    """
    _, size = vecs.shape
    P = np.eye(size)
    for vec in vecs:
        P = P - np.outer(vec, vec)
    return P


def get_hessian_projector(
    cart_coords,
    masses,
    cart_gradient=None,
    full=False,
    # These thresholds correspond to the 'gau' threshold
    max_grad_thresh=4.5e-4,
    rms_grad_thresh=3e-4,
):
    tr_vecs = get_trans_rot_vectors(cart_coords, masses=masses)
    vecs = tr_vecs

    # Only project out along gradient direction when we are NOT at a stationary point.
    if cart_gradient is not None:
        at_stationary_point = (np.abs(cart_gradient).max() <= max_grad_thresh) and (
            rms(cart_gradient) <= rms_grad_thresh
        )
        # See section II C in [1] or 2.2 in [2]
        if not at_stationary_point:
            P0 = get_projector(tr_vecs)
            # Project translation & rotation from gradient. This is usually already done
            # in the QC codes.
            mw_gradient = cart_gradient / np.sqrt(np.repeat(masses, 3))
            mw_gradient = P0 @ mw_gradient
            # Construct normalized vector along reaction coordinate ...
            rx_vec = mw_gradient / np.linalg.norm(mw_gradient)
            # ... and add it to the array that is used to construct the projector.
            vecs = np.concatenate(
                (
                    vecs,
                    rx_vec[None,],
                ),
                axis=0,
            )
    if full:
        P = get_projector(vecs)
    else:
        U, s, _ = np.linalg.svd(vecs.T)
        P = U[:, s.size :].T
    return P


def get_trans_rot_projector(cart_coords, masses, cart_gradient=None, full=False):
    warnings.warn(
        "'get_trans_rot_projector()' is deprecated. Please use 'get_hessian_projector() "
        "instead.",
        DeprecationWarning,
    )
    return get_hessian_projector(
        cart_coords, masses, cart_gradient=cart_gradient, full=full
    )
