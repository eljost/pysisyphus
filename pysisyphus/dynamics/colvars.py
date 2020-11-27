import abc
import types

import numpy as np
import autograd
import autograd.numpy as anp

from pysisyphus.intcoords.derivatives import dq_b, dq_a, dq_d


class Colvar(metaclass=abc.ABCMeta):
    def __init__(self, force_agrad=False):
        try:
            getattr(self, "_gradient")
        except AttributeError:
            force_agrad = True

        # Set autograd gradient method, if no _gradient forces, or if it
        # was not implemented.
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

    def _wilson_gradient(self, func, c3d):
        """Gradient of primitive internal w.r.t. Cartesians."""
        grad = np.zeros_like(c3d)
        grad[self.indices] = func(*c3d[self.indices].flatten()).reshape(-1, 3)
        return grad


class CVDistance(Colvar):
    def __init__(self, indices, **kwargs):
        self.indices = list(indices)
        self.i, self.j = self.indices
        super().__init__(**kwargs)

    def value(self, c3d):
        return anp.linalg.norm(c3d[self.i] - c3d[self.j])

    def _gradient(self, c3d):
        return self._wilson_gradient(dq_b, c3d)


class CVBend(Colvar):
    def __init__(self, indices, **kwargs):
        self.indices = list(indices)
        # Bonded like
        # i---j <- central atom
        #     |
        #     k
        self.i, self.j, self.k = self.indices
        super().__init__(**kwargs)

    def value(self, c3d):
        u_dash = c3d[self.i] - c3d[self.j]
        v_dash = c3d[self.k] - c3d[self.j]
        u_norm = anp.linalg.norm(u_dash)
        v_norm = anp.linalg.norm(v_dash)
        u = u_dash / u_norm
        v = v_dash / v_norm
        rad = anp.arccos(anp.dot(u, v))
        return rad

    def _gradient(self, c3d):
        return self._wilson_gradient(dq_a, c3d)


class CVDihedral(Colvar):
    def __init__(self, indices, **kwargs):
        self.indices = list(indices)
        # Bonded like
        # i---j
        #     |
        #     k---l
        self.i, self.j, self.k, self.l = self.indices
        super().__init__(**kwargs)

    def value(self, c3d):
        u_dash = c3d[self.i] - c3d[self.j]
        v_dash = c3d[self.k] - c3d[self.j]
        u_norm = anp.linalg.norm(u_dash)
        v_norm = anp.linalg.norm(v_dash)
        u = u_dash / u_norm
        v = v_dash / v_norm
        rad = anp.arccos(anp.dot(u, v))
        return rad

    def _gradient(self, c3d):
        return self._wilson_gradient(dq_d, c3d)
