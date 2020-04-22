import abc

import numpy as np


class Primitive(metaclass=abc.ABCMeta):

    def __init__(self, indices, periodic=False, calc_kwargs=None):
        self.indices = indices
        self.periodic = periodic
        if calc_kwargs is None:
            calc_kwargs = ()
        self.calc_kwargs = calc_kwargs

    @staticmethod
    def parallel(u, v, thresh=1e-4):
        dot = u.dot(v) / (np.linalg.norm(u) * np.linalg.norm(v))
        return (1 - abs(dot)) < thresh

    @abc.abstractmethod
    def _calculate(*, coords3d, indices, gradient, **kwargs):
        pass

    def calculate(self, coords3d, indices=None, gradient=False):
        if indices is None:
            indices = self.indices

        # Gather calc_kwargs
        calc_kwargs = {key: getattr(self, key) for key in self.calc_kwargs}

        return self._calculate(
                    coords3d=coords3d,
                    indices=indices,
                    gradient=gradient,
                    **calc_kwargs,
        )

    def __str__(self):
        return f"{self.__class__.__name__}({self.indices})"
