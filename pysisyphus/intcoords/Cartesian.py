import numpy as np

from pysisyphus.intcoords.Primitive import Primitive


class Cartesian(Primitive):

    def __init__(self, *args, **kwargs):
        kwargs["calc_kwargs"] = ("cart_axis",)
        super().__init__(*args, **kwargs)

    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        return 1

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, cart_axis=0):
        m,  = indices
        value = coords3d[m, cart_axis]
        if gradient:
            row = np.zeros_like(coords3d)
            row[m, cart_axis] = 1
            row = row.flatten()
            return value, row
        return value

    @staticmethod
    def _jacobian(coords3d, indices, cart_axis):
        return np.zeros(1)


class CartesianX(Cartesian):
    cart_axis = 0


class CartesianY(Cartesian):
    cart_axis = 1


class CartesianZ(Cartesian):
    cart_axis = 2
