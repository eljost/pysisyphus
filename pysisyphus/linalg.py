from math import cos, sin, sqrt
from typing import Callable, Literal

import numpy as np
import numpy.typing as npt


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
    coords: npt.NDArray[float],
    grad_func: Callable[[npt.NDArray[float]], npt.NDArray[float]],
    step_size: float = 1e-2,
    acc: Literal[2, 4] = 2,
) -> npt.NDArray[float]:
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
