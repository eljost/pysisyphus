import numpy as np


class Gaussian:
    def __init__(self, w=1, s=1, x0=None, cr_func=None):
        """
        See https://doi.org/10.1016/j.cpc.2018.02.017

        V(f(x)) = w * exp(-(f(x) - f0)**2 / (2*s**2))

        F(x) = -dV/dx = -dV/df * df/dx
        """
        # Gaussian height
        self._w = w
        # Standard deviation
        self._s = s
        # Center
        self.x0 = x0
        # Additional function for chain rule gradient
        if cr_func is None:
            cr_func = lambda x: 1.0
        self.cr_func = cr_func

        # Store some values, to avoid recalculation
        self.s2 = self.s ** 2
        self.one_over_2s2 = 1 / (2 * self.s2)
        self.minus_w_over_s2 = -self.w / self.s2

    @property
    def w(self):
        return self._w

    @property
    def s(self):
        return self._s

    def calculate(self, x, x0=None, gradient=False):
        """Return potential and gradient for Gaussian(s) centered at x0."""

        if x0 is None:
            x0 = self.x0
        diff = x - x0
        exp_ = np.exp(-(diff ** 2) * self.one_over_2s2)
        to_return = self.w * exp_.sum()
        if gradient:
            grad = self.minus_w_over_s2 * (diff * exp_ * self.cr_func(x)).sum(axis=0)
            # Append gradient
            to_return = to_return, grad
        return to_return

    def value(self, x, x0=None):
        return self.calculate(x, x0)

    def gradient(self, x, x0=None):
        _, grad = self.calculate(x, x0, gradient=True)
        return grad

    def eval(self, x, x0=None):
        return self.calculate(x, x0, gradient=True)
