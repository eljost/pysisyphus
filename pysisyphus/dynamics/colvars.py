import abc
import types

import numpy as np
import autograd
import autograd.numpy as anp

from pysisyphus.intcoords.derivatives import dq_b


class Colvar(metaclass=abc.ABCMeta):
    def __init__(self, force_agrad=False):
        try:
            _grad = getattr(self, "_gradient")
        except AttributeError:
            force_agrad = True

        if force_agrad:
            grad_func = autograd.grad(self.value)
            def wrapped(self, c3d):
                return grad_func(c3d)
            self._gradient = types.MethodType(wrapped, self)
        # Store a flag to indicate use of autograd
        self.agrad = force_agrad

    @abc.abstractmethod
    def value(self, c3d):
        pass

    def gradient(self, c3d):
        return self._gradient(c3d)

    def eval(self, c3d):
        return self.value(c3d), self.gradient(c3d)


class CVDistance(Colvar):
    def __init__(self, indices, **kwargs):
        self.indices = list(indices)
        self.i, self.j = self.indices
        super().__init__(**kwargs)

    def value(self, c3d):
        return anp.linalg.norm(c3d[self.i] - c3d[self.j])

    def _gradient(self, c3d):
        grad = np.zeros_like(c3d)
        grad[self.indices] = dq_b(*c3d[self.indices].flatten()).reshape(-1, 3)
        return grad
