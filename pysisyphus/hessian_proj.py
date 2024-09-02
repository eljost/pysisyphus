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


import dataclasses
from typing import Optional
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


@dataclasses.dataclass
class Orientation:
    """Orientation store, e.g., for a standard orientation.

    After instantiation the object performs a quick consistency check
    by trying to backtransform the re-oriented coordinates.

    Parameters
    ----------
    coords3d
        2d-array of shape (N, 3) holding re-oriented coordinates, e.g.,
        in standard orientation.
    com
        (Original) center-of-mass to restore original coordinates
        from coords3d.
    rot_mat
        2d-array of shape (3, 3) that was used to rotate the shifted
        original coordinates into the chosen orientation.
    coords3d_org
        2d-array of shape (N, 3) holding the original Cartesian coordinates.
    masses
        1d-arary of shape (N, ) holding atomic masses in amu.
    """

    coords3d: np.ndarray
    com: np.ndarray
    rot_mat: np.ndarray
    coords3d_org: np.ndarray
    masses: np.ndarray

    def __post_init__(self):
        # Quick sanity check
        tmp = np.einsum("ij,kj->ki", self.rot_mat.T, self.coords3d) + self.com
        np.testing.assert_allclose(tmp, self.coords3d_org, atol=1e-10)

    def rotate_gradient(self, gradient):
        gradient_rot = np.einsum("ij,kj->ki", self.rot_mat, gradient)
        return gradient_rot


def get_unit_orientation(coords3d, masses) -> Orientation:
    return Orientation(
        coords3d=coords3d.copy(),
        com=np.zeros(3),
        rot_mat=np.eye(3),
        coords3d_org=coords3d.copy(),
        masses=masses,
    )


def get_standard_orientation(coords3d, masses) -> Orientation:
    coords3d = coords3d.reshape(-1, 3)
    assert len(coords3d) == len(masses)
    # Center-of-mass
    com = 1 / masses.sum() * np.sum(coords3d * masses[:, None], axis=0)

    # Translate center of mass to origin
    coords3d_com = coords3d - com
    coords3d_aligned = coords3d_com.copy()
    rot_mat = np.eye(3)
    for _ in range(5):
        I = inertia_tensor(coords3d_aligned, masses)
        _, v = np.linalg.eigh(I)
        # Leave loop when principal axis are aligned
        if np.allclose(v, np.eye(3)):
            break
        # Otherwise align coordinates and updated rotation matrix
        coords3d_aligned = np.einsum("ij,kj->ki", v.T, coords3d_aligned)
        rot_mat = v.T @ rot_mat
    std_orient = Orientation(
        coords3d=coords3d_aligned,
        com=com,
        rot_mat=rot_mat,
        coords3d_org=coords3d.copy(),
        masses=masses,
    )
    return std_orient


def get_trans_rot_vectors(cart_coords, masses, rot_thresh=1e-5):
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
    # Drop vectors with vanishing norms/very small entries. Previous versions
    # tested the norm of these vectors. This was/is a bad idea, as the size of these
    # vectors grows with the number of atoms. Now we use the rms value instead.
    rot_vec_rms = np.sqrt(np.mean(rot_vecs**2, axis=1))
    mask = rot_vec_rms > rot_thresh
    rot_vecs = rot_vecs[mask]
    # Vectors are in rows after concatenate call
    tr_vecs = np.concatenate((trans_vecs, rot_vecs), axis=0)
    # Orthonormalize vectors
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
    cart_coords: np.ndarray,
    masses: np.ndarray,
    cart_gradient: Optional[np.ndarray] = None,
    full: bool = False,
    # These thresholds correspond to the 'gau' threshold
    max_grad_thresh: float = 4.5e-4,
    rms_grad_thresh: float = 3e-4,
    use_std_orient: bool = True,
):
    """Get matrix to project translation and rotation from the Hessian.

    Parameters
    ----------
    cart_coords
        1d array of shape (3N, ) holding Cartesian coordinates.
    masses
        1d arary of shape (N, ) holding atomic masses.
    cart_gradient
        Optional 1d array of shape (3N, ) containing a Cartesian gradient
        that will also be projected out of the Hessian.
    full
        Boolean flag that controls the shape of the projector. If true,
        the projector will be of shape (3N, 3N), otherwise it will be of
        shape (3N - 6 (5), 3N - 6 (5)).
    max_grad_thresh
        Positive floating point number. Criterion used to determine if we are
        at a stationary point. If we are at a statioanry point, the gradient
        direction won't be included in the projector.
    rms_grad_thresh
        See max_grad_thresh.
    use_std_orient
        Boolean flag that controls the use of a (temporary) standard orientation.
        If set to False linear molecules can lack a vibration.
    """
    orient_func = get_standard_orientation if use_std_orient else get_unit_orientation
    orient = orient_func(cart_coords.reshape(-1, 3), masses)
    # Continue with (possibly) re-oriented coordinates
    cart_coords = orient.coords3d.flatten()
    if cart_gradient is not None:
        cart_gradient = orient.rotate_gradient(cart_gradient)

    # Calculate vectors to project out translation and rotation
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
    # Bring vectors back into original orientation. Note the transposed rotation matrix,
    # to restore the original orientation.
    natoms = len(masses)
    eye = np.eye(natoms)
    rot_mat = np.kron(eye, orient.rot_mat.T)
    P = (rot_mat @ P.T).T
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
