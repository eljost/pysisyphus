import abc

import numpy as np


class Primitive(metaclass=abc.ABCMeta):

    def __init__(self, indices, periodic=False):
        self.indices = indices
        self.periodic = periodic

    @staticmethod
    def parallel(u, v, thresh=1e-4):
        dot = u.dot(v) / (np.linalg.norm(u) * np.linalg.norm(v))
        return (1 - abs(dot)) < thresh

    @abc.abstractmethod
    def _calculate(*, coords3d, indices, gradient):
        pass

    def calculate(self, coords3d, indices=None, gradient=False):
        if indices is None:
            indices = self.indices

        return self._calculate(
                    coords3d=coords3d,
                    indices=indices,
                    gradient=gradient
        )
