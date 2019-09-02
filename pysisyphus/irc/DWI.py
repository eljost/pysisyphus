#!/usr/bin/env python3

# [1] https://aip.scitation.org/doi/10.1063/1.1724823
#     Hratchian, 2004
# [2] https://aip.scitation.org/doi/pdf/10.1063/1.475419?class=pdf
#     Thompson, 1998

from collections import deque

import numpy as np


def taylor(energy, gradient, hessian, step):
    return energy + step @ gradient + 0.5 * step @ hessian @ step


class DWI:

    def __init__(self, n=4, maxlen=2):
        self.n = int(n)
        assert self.n > 0
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

    def interpolate(self, at_coords):
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

        return E_dwi
