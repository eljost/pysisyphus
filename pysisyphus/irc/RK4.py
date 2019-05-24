#!/usr/bin/env python3

import numpy as np

from pysisyphus.irc.IRC import IRC


# [1] https://aip.scitation.org/doi/pdf/10.1063/1.462674?class=pdf
# [2] https://www.sciencedirect.com/science/article/pii/S0009261406015673


class RK4(IRC):

    def get_k(self, mw_coords):
        self.mw_coords = mw_coords
        grad = self.mw_gradient
        direction = -grad / np.linalg.norm(grad)
        k = self.step_length * direction
        return k

    def step(self):
        # Eq. (25 - 29) in [1]
        mw_coords_1 = self.mw_coords.copy()
        k1 = self.get_k(mw_coords_1)

        mw_coords_2 = mw_coords_1 + 0.5*k1
        k2 = self.get_k(mw_coords_2)

        mw_coords_3 = mw_coords_1 + 0.5*k2
        k3 = self.get_k(mw_coords_3)

        mw_coords_4 = mw_coords_1 + k3
        k4 = self.get_k(mw_coords_4)

        step = (k1 + 2*k2 + 2*k3 + k4) / 6
        step_norm = np.linalg.norm(step)
        self.log(f"norm(step)={step_norm:.6f}")
        self.mw_coords = mw_coords_1 + step
