# [1] http://dx.doi.org/10.1063/1.4952956
#     Lee-Ping Wang, 2016

import numpy as np

from pysisyphus.intcoords.Primitive import Primitive


class Translation(Primitive):
    """See (II. Theory) in [1], Eq. (2)"""

    def __init__(self, *args, **kwargs):
        kwargs["calc_kwargs"] = ("cart_axis",)
        super().__init__(*args, **kwargs)

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        return 1

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, cart_axis=0):
        value = coords3d[indices, cart_axis].mean()
        if gradient:
            row = np.zeros_like(coords3d)
            row[indices, cart_axis] = 1 / len(indices)
            row = row.flatten()
            return value, row
        return value

    @staticmethod
    def _jacobian(coords3d, indices, cart_axis):
        size = len(indices) * 3
        return np.zeros(size * size)


class TranslationX(Translation):
    cart_axis = 0


class TranslationY(Translation):
    cart_axis = 1


class TranslationZ(Translation):
    cart_axis = 2
