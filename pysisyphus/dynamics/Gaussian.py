import numpy as np


class DummyColvar:
    def eval(self, x):
        assert np.isscalar(
            x
        ), f"DummyColvar can only be used for scalar inputs, but got {x}!"
        return x, 1


class Gaussian:
    def __init__(self, w=1, s=1, x0=None, colvar=None):
        """
        See https://doi.org/10.1016/j.cpc.2018.02.017

        V(f(x)) = w * exp(-(f(x) - f0)**2 / (2*s**2))

        F(x) = -dV/dx = -dV/df * df/dx
        """
        # Gaussian height w and standard deviation s are read-only properties.
        # This way they can't be accidentally altered, which would invalidate the
        # values precomputed below.
        self._w = w
        self._s = s
        self.x0 = np.ravel(x0)

        # Collective variable
        if colvar is None:
            colvar = DummyColvar()
        self.colvar = colvar

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

    def calculate(self, coords, x0=None, gradient=False):
        """Return potential and gradient for Gaussian(s) centered at x0."""

        if x0 is None:
            x0 = self.x0
        x, cr_grad = self.colvar.eval(coords)
        diff = x - np.atleast_1d(x0)
        exp_ = np.exp(-(diff ** 2) * self.one_over_2s2)
        to_return = self.w * exp_.sum()

        if gradient:
            try:
                grad = np.einsum("i,i,ijk->jk", diff, exp_, cr_grad[None, :, :])
            except TypeError:
                grad = diff * exp_ * cr_grad

            # Finalize & append gradient
            grad *= self.minus_w_over_s2
            to_return = to_return, grad
        return to_return

    def value(self, coords, x0=None):
        return self.calculate(coords, x0)

    def gradient(self, coords, x0=None):
        _, grad = self.calculate(coords, x0, gradient=True)
        return grad

    def eval(self, coords, x0=None):
        return self.calculate(coords, x0, gradient=True)
