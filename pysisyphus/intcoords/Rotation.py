# [1] http://dx.doi.org/10.1063/1.4952956
#     Lee-Ping Wang, 2016

import numpy as np

from pysisyphus.intcoords.Primitive import Primitive
from pysisyphus.linalg import eigvec_grad, rot_quaternion


def compare_to_geometric(c3d, ref_c3d, dR, dF, dqdx, dvdx, atol=1e-14):
    from geometric.rotate import get_R_der, get_F_der, get_q_der, get_expmap_der

    dR_ref = get_R_der(c3d, ref_c3d)
    np.testing.assert_allclose(dR, dR_ref)
    dF_ref = get_F_der(c3d, ref_c3d)
    np.testing.assert_allclose(dF.reshape(-1, 3, 4, 4), dF_ref)
    dq_ref = get_q_der(c3d, ref_c3d)
    np.testing.assert_allclose(dqdx.T.flatten(), dq_ref.flatten(), atol=atol)
    dvdx_ref = get_expmap_der(c3d, ref_c3d)
    np.testing.assert_allclose(dvdx.T.flatten(), dvdx_ref.flatten(), atol=atol)


class Rotation(Primitive):
    """See (II. Theory) in [1], Eq. (3) - (14)"""

    index = None

    def __init__(self, indices, *args, ref_coords3d, **kwargs):
        kwargs["cache"] = False
        kwargs["calc_kwargs"] = ("index", "ref_coords3d")
        super().__init__(indices, *args, **kwargs)

        self.ref_coords3d = ref_coords3d.reshape(-1, 3).copy()

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        return 1

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, index=0, ref_coords3d=None):
        w, v_, c3d, ref_c3d = rot_quaternion(coords3d[indices], ref_coords3d[indices])

        quat = v_[:, -1]
        # Eigenvector sign is ambigous. Force first item to be positive,
        # similar to geomeTRIC code.
        if quat[0] < 0.0:
            quat *= -1

        # Eq. (8) in [1].
        # v = 2 * q_i * (cos⁻¹(q_0) / sqrt(1 - q_0 ** 2)
        #
        # As q_0 approaches 1, the denominator becomes very small, and dividing
        # by this small number results in numerical instability.
        #
        # According to wolframalpha v(q_0) limit approaches 2 for q_0 = 1.
        #
        #   input: limit of (2 * arccos(x) / sqrt(1-x**2))
        #   output: lim v(x) for x -> 1 becomes 2.
        q0 = quat[0]
        if abs(q0 - 1.0) <= 1e-8:
            prefac = 2 - 2 / 3 * (q0 - 1)
            dvdq0 = -2 / 3
        else:
            arccos_q0 = np.arccos(q0)
            diff = 1 - q0 ** 2
            prefac = 2 * arccos_q0 / np.sqrt(diff)
            dvdq0 = 2 * q0 * arccos_q0 / diff ** 1.5 - 2 / diff

        # Exponential map
        v = prefac * quat[1:]

        if gradient:
            # Gradient of correlation matrix
            y1, y2, y3 = ref_c3d.T
            dR = np.zeros((*c3d.shape, 3, 3))
            dR[:, 0, 0, 0] = y1
            dR[:, 0, 0, 1] = y2
            dR[:, 0, 0, 2] = y3
            #
            dR[:, 1, 1, 0] = y1
            dR[:, 1, 1, 1] = y2
            dR[:, 1, 1, 2] = y3
            #
            dR[:, 2, 2, 0] = y1
            dR[:, 2, 2, 1] = y2
            dR[:, 2, 2, 2] = y3
            dR11, dR12, dR13, dR21, dR22, dR23, dR31, dR32, dR33 = dR.reshape(-1, 9).T

            # Gradient of F matrix. Construct full matrix, as we have to do a dot
            # product later on.
            dF = np.zeros((ref_c3d.size, 4, 4))
            dF[:, 0, 0] = dR11 + dR22 + dR33
            dF[:, 0, 1] = dR23 - dR32
            dF[:, 0, 2] = dR31 - dR13
            dF[:, 0, 3] = dR12 - dR21
            #
            dF[:, 1, 0] = dF[:, 0, 1]
            dF[:, 1, 1] = dR11 - dR22 - dR33
            dF[:, 1, 2] = dR12 + dR21
            dF[:, 1, 3] = dR13 + dR31
            #
            dF[:, 2, 0] = dF[:, 0, 2]
            dF[:, 2, 1] = dF[:, 1, 2]
            dF[:, 2, 2] = -dR11 + dR22 - dR33
            dF[:, 2, 3] = dR23 + dR32
            #
            dF[:, 3, 0] = dF[:, 0, 3]
            dF[:, 3, 1] = dF[:, 1, 3]
            dF[:, 3, 2] = dF[:, 2, 3]
            dF[:, 3, 3] = -dR11 - dR22 + dR33

            # Quaternion gradient
            dqdx = eigvec_grad(w, v_, ind=-1, mat_grad=dF)

            dvdq = np.zeros((3, 4))
            dvdq[:, 0] = dvdq0 * quat[1:]
            dvdq[:, 1:] = np.diag((prefac, prefac, prefac))

            # Gradient of exponential map from chain rule.
            # See bottom-left on 214108-3 in [1], after Eq. (11).
            dvdx = np.einsum("ji,ik->jk", dvdq, dqdx)

            # compare_to_geometric(c3d, ref_c3d, dR, dF, dqdx, dvdx)
            row = np.zeros_like(coords3d)
            if index is None:
                return v, dvdx.reshape(3, -1)
            row[indices] = dvdx[index].reshape(-1, 3)
            return v[index], row.flatten()
        return v[index]

    @staticmethod
    def _jacobian(coords3d, indices, index=0, ref_coords3d=None):
        """Not implemented!"""
        size = len(indices) * 3
        return np.zeros(size * size)


class RotationA(Rotation):
    index = 0


class RotationB(Rotation):
    index = 1


class RotationC(Rotation):
    index = 2
