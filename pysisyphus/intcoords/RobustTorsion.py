# [1] Transition-State Optimization Methods using Internal Coordinates
#     https://macsphere.mcmaster.ca/handle/11375/15450?mode=full
#     Rabi, 2014, Phd Thesis

import numpy as np

from pysisyphus.intcoords.Primitive import Primitive
from pysisyphus.intcoords.derivatives import (
    q_rd1,
    dq_rd1,
    d2q_rd1,
    q_rd2,
    dq_rd2,
    d2q_rd2,
)


class RobustTorsion1(Primitive):
    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        return 1.0

    @staticmethod
    def _calculate(coords3d, indices, gradient=False):
        args = coords3d[indices].flatten()
        val = q_rd1(*args)

        if gradient:
            row = np.zeros_like(coords3d)
            grad = dq_rd1(*args).reshape(-1, 3)
            row[indices] = grad
            return val, row.flatten()
        return val

    @staticmethod
    def _jacobian(coords3d, indices):
        return d2q_rd1(*coords3d[indices].flatten())


class RobustTorsion2(Primitive):
    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        return 1.0

    @staticmethod
    def _calculate(coords3d, indices, gradient=False):
        args = coords3d[indices].flatten()
        val = q_rd2(*args)

        if gradient:
            row = np.zeros_like(coords3d)
            grad = dq_rd2(*args).reshape(-1, 3)
            row[indices] = grad
            return val, row.flatten()
        return val

    @staticmethod
    def _jacobian(coords3d, indices):
        return d2q_rd2(*coords3d[indices].flatten())
