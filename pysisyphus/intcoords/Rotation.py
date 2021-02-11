# [1] http://dx.doi.org/10.1063/1.4952956
#     Lee-Ping Wang, 2016

import numpy as np

from pysisyphus.intcoords.Primitive import Primitive


class Rotation(Primitive):
    """See (II. Theory) in [1], Eq. (3) - (14)"""

    index = None

    def __init__(self, *args, ref_coords3d, **kwargs):
        super().__init__(*args, **kwargs)

        self.calc_kwargs = ("index", "ref_coords3d")
        self.ref_coords3d = ref_coords3d.reshape(-1, 3).copy()

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        return 1

    @staticmethod
    def to_origin(coords3d, indices):
        return coords3d[indices] - coords3d[indices].mean(axis=0)

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, index=0, ref_coords3d=None):
        # Translate to origin by removing centroid
        c3d = Rotation.to_origin(coords3d, indices)
        ref_c3d = Rotation.to_origin(ref_coords3d, indices)

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
        w, v = np.linalg.eigh(F, UPLO="U")
        # Quaternion corresponds to biggest (last) eigenvalue.
        # np.linalg.eigh already returns sorted eigenvalues.
        quat = v[:, -1]
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
            prefac = 2
        else:
            prefac = 2 * np.arccos(q0) / np.sqrt(1 - q0 ** 2)
        # Exponential map
        v = prefac * quat[1:]

        if gradient:
            """See https://github.com/google/jax/issues/2748"""
            raise Exception("Not implemented!")

        return v[index]

    @staticmethod
    def _jacobian(coords3d, indices):
        raise Exception("Not implemented!")


class RotationA(Rotation):
    index = 0


class RotationB(Rotation):
    index = 1


class RotationC(Rotation):
    index = 2
