# See
# https://www.jwwalker.com/pages/angle-between-vectors.html

import numpy as np

from pysisyphus.intcoords.Bend import Bend
from pysisyphus.intcoords.derivatives import q_a2, dq_a2, d2q_a2


class Bend2(Bend):

    @staticmethod
    def _calculate(coords3d, indices, gradient=False):
        args = coords3d[indices].flatten()
        val = q_a2(*args)
        if gradient:
            grad = dq_a2(*args).reshape(-1, 3)
            row = np.zeros_like(coords3d)
            row[indices] = grad
            return val, row.flatten()
        return  val

    @staticmethod
    def _jacobian(coords3d, indices):
        return d2q_a2(*coords3d[indices].flatten())
