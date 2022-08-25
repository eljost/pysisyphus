from math import cos, sin, sqrt
from typing import Callable, Literal, Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.spatial.transform import Rotation
from scipy.linalg.lapack import dpstrf


def gram_schmidt(vecs, thresh=1e-8):
    def proj(v1, v2):
        return v1.dot(v2) / v1.dot(v1)

    ortho = [
        vecs[0],
    ]
    for v1 in vecs[1:]:
        tmp = v1.copy()
        for v2 in ortho:
            tmp -= proj(v2, v1) * v2
        norm = np.linalg.norm(tmp)
        # Don't append linear dependent vectors, as their norm will be
        # near zero. Renormalizing them to unity would lead to numerical
        # garbage and to erronous results later on, when we orthgonalize
        # against this 'arbitrary' vector.
        if norm <= thresh:
            continue
        ortho.append(tmp / norm)
    return np.array(ortho)


def orthogonalize_against(mat, vecs, max_cycles=5, thresh=1e-10):
    """Orthogonalize rows of 'mat' against rows in 'vecs'

    Returns a (modified) copy of mat.
    """

    omat = mat.copy()
    for _ in range(max_cycles):
        max_overlap = 0.0
        for row in omat:
            for ovec in vecs:
                row -= ovec.dot(row) * ovec

            # Remaining overlap
            overlap = np.sum(np.dot(vecs, row))
            max_overlap = max(overlap, max_overlap)

        if max_overlap < thresh:
            break
    else:
        raise Exception(f"Orthogonalization did not succeed in {max_cycles} cycles!")

    return omat


def perp_comp(vec, along):
    """Return the perpendicular component of vec along along."""
    return vec - vec.dot(along) * along


def make_unit_vec(vec1, vec2):
    """Return unit vector pointing from vec2 to vec1."""
    diff = vec1 - vec2
    return diff / np.linalg.norm(diff)


def svd_inv(array, thresh, hermitian=False):
    U, S, Vt = np.linalg.svd(array, hermitian=hermitian)
    keep = S > thresh
    S_inv = np.zeros_like(S)
    S_inv[keep] = 1 / S[keep]
    return Vt.T.dot(np.diag(S_inv)).dot(U.T)


def get_rot_mat(abc=None):
    # Euler angles
    if abc is None:
        abc = np.random.rand(3) * np.pi * 2
    a, b, c = abc
    R = np.array(
        (
            (
                cos(a) * cos(b) * cos(c) - sin(a) * sin(c),
                -cos(a) * cos(b) * sin(c) - sin(a) * cos(c),
                cos(a) * sin(b),
            ),
            (
                sin(a) * cos(b) * cos(c) + cos(a) * sin(c),
                -sin(a) * cos(b) * sin(c) + cos(a) * cos(c),
                sin(a) * sin(b),
            ),
            (-sin(b) * cos(c), sin(b) * sin(c), cos(b)),
        )
    )
    return R


def get_rot_mat_for_coords(coords3d_1, coords3d_2):
    coords3d_1 = coords3d_1.copy()
    coords3d_2 = coords3d_2.copy()
    centroid_1 = coords3d_1.mean(axis=0)
    centroid_2 = coords3d_2.mean(axis=0)
    coords3d_1 -= centroid_1[None, :]
    coords3d_2 -= centroid_2[None, :]

    tmp = coords3d_2.T.dot(coords3d_1)
    U, W, Vt = np.linalg.svd(tmp)
    rot_mat = U.dot(Vt)
    if np.linalg.det(rot_mat) < 0:
        U[:, -1] *= -1
        rot_mat = U.dot(Vt)
    return coords3d_1, coords3d_2, rot_mat


def eigvec_grad(w, v, ind, mat_grad):
    """Gradient of 'ind'-th eigenvector.

    dv_i / dx_i = (w_i*I - mat)⁻¹ dmat/dx_i v_i
    """
    eigval = w[ind]
    eigvec = v[:, ind]

    w_diff = eigval - w
    w_inv = np.divide(
        1.0, w_diff, out=np.zeros_like(w_diff).astype(float), where=w_diff != 0.0
    )
    assert np.isfinite(w_inv).all()
    pinv = v.dot(np.diag(w_inv)).dot(v.T)
    pinv.dot(mat_grad)

    wh = pinv.dot(mat_grad).dot(eigvec)
    return wh


def cross3(a, b):
    """10x as fast as np.cross for two 1d arrays of size 3."""
    return np.array(
        (
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        )
    )


def norm3(a):
    """5x as fas as np.linalg.norm for a 1d array of size 3."""
    return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])


def finite_difference_hessian(
    coords: NDArray[float],
    grad_func: Callable[[NDArray[float]], NDArray[float]],
    step_size: float = 1e-2,
    acc: Literal[2, 4] = 2,
) -> NDArray[float]:
    """Numerical Hessian from central finite gradient differences.

    See central differences in
      https://en.wikipedia.org/wiki/Finite_difference_coefficient
    for the different accuracies.
    """
    accuracies = {
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
        step[i] = step_size

        def get_grad(factor, displ):
            displ_coords = coords + step * displ
            grad = grad_func(displ_coords)
            return factor * grad

        grads = [get_grad(factor, displ) for factor, displ in coeffs]
        fd = np.sum(grads, axis=0) / step_size
        fd_hessian[i] = fd

    # Symmetrize
    fd_hessian = (fd_hessian + fd_hessian.T) / 2
    return fd_hessian


def rot_quaternion(coords3d, ref_coords3d):
    # Translate to origin by removing centroid
    c3d = coords3d - coords3d.mean(axis=0)
    ref_c3d = ref_coords3d - ref_coords3d.mean(axis=0)

    # Setup correlation matrix
    R = c3d.T.dot(ref_c3d)

    # Setup F matrix, Eq. (6) in [1]
    F = np.zeros((4, 4))
    R11, R12, R13, R21, R22, R23, R31, R32, R33 = R.flatten()
    # Fill only upper triangular part.
    F[0, 0] = R11 + R22 + R33
    F[0, 1] = R23 - R32
    F[0, 2] = R31 - R13
    F[0, 3] = R12 - R21
    #
    F[1, 1] = R11 - R22 - R33
    F[1, 2] = R12 + R21
    F[1, 3] = R13 + R31
    #
    F[2, 2] = -R11 + R22 - R33
    F[2, 3] = R23 + R32
    #
    F[3, 3] = -R11 - R22 + R33

    # Eigenvalues, eigenvectors of upper triangular part.
    w, v_ = np.linalg.eigh(F, UPLO="U")

    # Quaternion corresponds to biggest (last) eigenvalue.
    # np.linalg.eigh already returns sorted eigenvalues.
    return w, v_, c3d, ref_c3d


def quaternion_to_rot_mat(q):
    q0, q1, q2, q3 = q
    q_ = q0 ** 2 - (q1 ** 2 + q2 ** 2 + q3 ** 2)
    R = np.zeros((3, 3))
    R[0, 0] = q_ + 2 * q1 ** 2
    R[0, 1] = 2 * (q1 * q2 - q0 * q3)
    R[0, 2] = 2 * (q1 * q3 - q0 * q2)
    R[1, 0] = 2 * (q1 * q2 + q0 * q3)
    R[1, 1] = q_ + 2 * q2 ** 2
    R[1, 2] = 2 * (q2 * q3 - q0 * q1)
    R[2, 0] = 2 * (q1 * q3 - q0 * q2)
    R[2, 1] = 2 * (q2 * q3 + q0 * q1)
    R[2, 2] = q_ + 2 * q3 ** 2
    return R


def rmsd_grad(
    coords3d: NDArray[float], ref_coords3d: NDArray[float], offset: float = 1e-9
) -> Tuple[float, NDArray[float]]:
    """RMSD and gradient between two sets of coordinates from quaternions.

    The gradient is given w.r.t. the coordinates of 'coords3d'.

    Python adaption of
        ls_rmsd.f90
    from the xtb repository of the Grimme group, which in turn implements
        [1] https://doi.org/10.1002/jcc.20110
    .

    Parameters
    ----------
    coords3d
        Coordinate array of shape (N, 3) with N denoting the number of atoms.
    ref_coords3d
        Reference coordinates.
    offset
        Small floating-point number that is added to the RMSD, to avoid division
        by zero.

    Returns
    -------
    rmsd
        RMSD value.
    rmsd_grad
        Gradient of the RMSD value w.r.t. to 'coords3d'.
    """
    assert coords3d.shape == ref_coords3d.shape

    w, v_, c3d, ref_c3d = rot_quaternion(coords3d, ref_coords3d)
    quat = v_[:, -1]
    eigval = w[-1]

    atom_num = coords3d.shape[0]
    x_norm = np.linalg.norm(c3d) ** 2
    y_norm = np.linalg.norm(ref_c3d) ** 2
    rmsd = np.sqrt(max(0.0, ((x_norm + y_norm) - 2.0 * eigval)) / (atom_num)) + offset

    # scipy expects the quaternion
    #     a + b*i + c*j + d*k
    # as (i, j, k, a), that is the scalar component must come last.
    # Currently we have (a, i, j, k), so we have to rearrange before passing the
    # quaternion to scipy.
    rot = Rotation.from_quat((*quat[1:], quat[0]))
    U = rot.as_matrix()

    grad = (c3d - ref_c3d @ U) / (rmsd * atom_num)
    return rmsd, grad


def pivoted_cholesky(A: NDArray, tol: float = -1.0):
    """Cholesky factorization a real symmetric positive semidefinite matrix.
    Cholesky factorization is carried out with full pivoting.

    Adapated from PySCF.

    P.T * A * P = L * L.T

    Parameters
    ----------
    A
        Matrix to be factorized.
    tol
        User defined tolerance, as outlined in the LAPACK docs.

    Returns
    -------
    L
        Lower or upper triangular matrix.
    piv
        Pivot vectors, starting at 0.
    rank
        Rank of the factoirzed matrix.
    """

    N = A.shape[0]
    L, piv, rank, info = dpstrf(A, tol=tol, lower=True)
    piv -= 1  # LAPACK returns 1-based indices
    assert info >= 0
    L[np.triu_indices(N, k=1)] = 0
    L[:, rank:] = 0
    return L, piv, rank
