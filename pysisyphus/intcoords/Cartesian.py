import numpy as np

from pysisyphus.intcoords.Primitive import Primitive


class Cartesian(Primitive):
    @staticmethod
    def _weight(atoms, coords3d, indices, f_damping):
        return 1

    @staticmethod
    def _calculate(coords3d, indices, gradient=False, cart_axis=0):
        value = coords3d[indices, cart_axis]
        if gradient:
            row = np.zeros_like(coords3d)
            row[indices, cart_axis] = 1
            row = row.flatten()
            return value, row
        return value

    @staticmethod
    def _jacobian(coords3d, indices):
        return np.zeros(1)


class CartesianX(Cartesian):
    cart_axis = 0


class CartesianY(Cartesian):
    cart_axis = 1


class CartesianZ(Cartesian):
    cart_axis = 2
