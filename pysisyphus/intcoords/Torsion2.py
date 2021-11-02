# See
# https://en.wikipedia.org/wiki/Dihedral_angle#In_polymer_physics

import numpy as np

from pysisyphus.intcoords import Torsion
from pysisyphus.intcoords.derivatives import q_d2, dq_d2, d2q_d2


class Torsion2(Torsion):
    @staticmethod
    def _calculate(coords3d, indices, gradient=False):
        args = coords3d[indices].flatten()
        val = q_d2(*args)

        if gradient:
            row = np.zeros_like(coords3d)
            grad = dq_d2(*args).reshape(-1, 3)
            row[indices] = grad
            return val, row.flatten()
        return val

    @staticmethod
    def _jacobian(coords3d, indices):
        return d2q_d2(*coords3d[indices].flatten())
