#!/usr/bin/env python3

# [1] https://aip.scitation.org/doi/10.1063/1.1724823
#     Hratchian, 2004
# [2] https://aip.scitation.org/doi/pdf/10.1063/1.475419?class=pdf
#     Thompson, 1998

from collections import deque

import numpy as np


def taylor(energy, gradient, hessian, step):
    """Taylor series expansion of the energy to second order."""
    return energy + step @ gradient + 0.5 * step @ hessian @ step


def taylor_grad(gradient, hessian, step):
    """Gradient of a Taylor series expansion of the energy to second order."""
    return gradient + hessian @ step


class DWI:

    def __init__(self, n=4, maxlen=2):
        """Distance weighted interpolation."""

        self.n = int(n)
        assert self.n > 0
        assert (self.n % 2) == 0
        self.maxlen = maxlen
        assert self.maxlen == 2, \
            "Right now only maxlen=2 is supported!"

        # Using FIFO deques for easy updating of the lists
        self.coords = deque(maxlen=self.maxlen)
        self.energies = deque(maxlen=self.maxlen)
        self.gradients = deque(maxlen=self.maxlen)
        self.hessians = deque(maxlen=self.maxlen)
    
    def update(self, coords, energy, gradient, hessian):
        self.coords.append(coords)
        self.energies.append(energy)
        self.gradients.append(gradient)
        self.hessians.append(hessian)

        assert len(self.coords) == len(self.energies) \
               == len(self.gradients) == len(self.hessians)

    def interpolate(self, at_coords, gradient=False):
        """See [1] Eq. (25) - (29)"""
        c1, c2 = self.coords

        dx1 = at_coords - c1
        dx2 = at_coords - c2

        dx1_norm = np.linalg.norm(dx1)
        dx2_norm = np.linalg.norm(dx2)
        dx1_norm_n = dx1_norm**self.n
        dx2_norm_n = dx2_norm**self.n

        denom = dx1_norm**self.n + dx2_norm**self.n
        w1 = dx2_norm_n / denom
        w2 = dx1_norm_n / denom

        e1, e2 = self.energies
        g1, g2 = self.gradients
        h1, h2 = self.hessians

        t1 = taylor(e1, g1, h1, dx1)
        t2 = taylor(e2, g2, h2, dx2)

        E_dwi = w1*t1 + w2*t2

        if not gradient:
            return E_dwi

        t1_grad = taylor_grad(g1, h1, dx1)
        t2_grad = taylor_grad(g2, h2, dx2)

        # The gradient of dx2_norm_n w.r.t the coordinates is formulated with
        # **2n instead of **n, so the square root can be easily reduced.
        # sqrt(x)**2n = x**(1/2)**2n = x**n
        #
        # Thats why we have to do following calculations with the half the value
        # of n.
        n_2 = self.n // 2
        dx1_norm_n_grad = 2 * n_2 * dx1_norm**(2*n_2-2) * dx1
        dx2_norm_n_grad = 2 * n_2 * dx2_norm**(2*n_2-2) * dx2
        w1_grad = (dx2_norm_n_grad*dx1_norm_n  - dx1_norm_n_grad*dx2_norm_n) / denom**2
        w2_grad = -w1_grad

        # E_dwi = w1(x)*T1(x) + w2(x)*T2(x)
        #
        # dE_DWI / dx = dw1(x)*T1(x) + w1(x)*dT1(x) + dw2(x)*T2(x) + w2*dT2(x)
        grad_dwi = w1_grad*t1 + w1*t1_grad + w2_grad*t2 + w2*t2_grad

        return E_dwi, grad_dwi
