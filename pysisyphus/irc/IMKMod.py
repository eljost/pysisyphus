#!/usr/bin/env python3

import warnings

import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.irc.IRC import IRC


# [1] http://pubs.acs.org/doi/pdf/10.1021/ja00295a002
# [2] https://math.stackexchange.com/questions/1882772/curve-fitting-with-derivatives


class IMKMod(IRC):

    def __init__(self, geometry, line_step_size=0.025, **kwargs):
        super().__init__(geometry, **kwargs)

        self.line_step_size = line_step_size
        self.line_step_thresh = 2*self.line_step_size

    def fit_parabola(self, x, y):
        y_str = np.array2string(np.array(y), precision=6)
        self.log(f"Fitting parabola with x={x} and y={y_str}")
        fit = np.polyfit(x, y, deg=2)
        parabola = np.poly1d(fit)
        minimum = parabola.deriv().r
        # We are only looking for a real minimum
        real_minimum = (minimum[minimum.imag==0].real)[0]
        self.log(f"Found minimum of fitted parabola at x={real_minimum:.4f}")

        return real_minimum

    def step(self):
        # Get gradient at Q0
        mw_coords_0 = self.mw_coords.copy()
        grad_0 = self.mw_gradient
        energy_0 = self.energy

        grad_0_norm = np.linalg.norm(grad_0)
        # Step direction is against the gradient to the pivot point Q1.
        step = -self.step_length * grad_0 / grad_0_norm
        mw_coords_1 = self.mw_coords + step
        self.mw_coords = mw_coords_1

        energy_1 = self.energy
        if energy_1 > energy_0:
            self.log("Initial displacement lead to higher energy. "
                     "Searching for better step with a shorter displacement."
            )
            corr_step = 0.5*step
            mw_coords_1_corr = mw_coords_0 + corr_step
            self.mw_coords = mw_coords_1_corr
            energy_1_corr = self.energy
            x_corr = np.array((0, 0.5, 1)) * self.step_length
            y_corr = (energy_0, energy_1_corr, energy_1)
            min_ = self.fit_parabola(x_corr, y_corr)
            step_corr = min_ * -grad_0 / grad_0_norm
            mw_coords_1 = mw_coords_0 + step_corr
            self.mw_coords = mw_coords_1

        # Get gradient at Q1
        grad_1 = self.mw_gradient
        energy_1 = self.energy
        grad_1_norm = np.linalg.norm(grad_1)

        # Determine bisector d
        d = -grad_0/grad_0_norm + grad_1/grad_1_norm
        D = d / np.linalg.norm(d)

        # Take a step along D
        step_2 = self.line_step_size * D
        mw_coords_2 = mw_coords_1 + step_2
        self.mw_coords = mw_coords_2
        energy_2 = self.energy
        if energy_2 > energy_1:
            factor = 0.5
        else:
            factor = 2
        step_3 = self.line_step_size * factor * D
        mw_coords_3 = mw_coords_1 + step_3
        self.mw_coords = mw_coords_3
        energy_3 = self.energy
        # x = (0, 1, factor)
        x = np.array((0, 1, factor)) * self.line_step_size
        y = (energy_1, energy_2, energy_3)
        min_ = self.fit_parabola(x, y)
        step_corr = min_ * D
        mw_coords_4 = mw_coords_1 + step_corr
        self.mw_coords = mw_coords_4
        return

        # # Project g1 along D
        # grad_1_proj = grad_1.dot(D)
        # """Fit parabel with 2 energies and 1 projected gradient.
        # See [2] on how to do this. We want to determine the coefficients
        # for the polynom f(x) = a0 + a1*x + a2*x**2. This can be done by
        # solving A*b = y with b being the coefficient vector.
        # The system matrix A is given by:
            # A = (
                # (1, x11, x11**2),
                # (1, x12, x12**2),
                # (0, 1, 2*x21)
            # )
        # energy_1 and and grad_1 are evaluated at the same point x11 = 0,
        # so only x12 is != 0."""
        # A = (
	    # (1, 0, 0),
	    # (1, self.line_step_size, self.line_step_size**2),
	    # (0, 1, 0)
	# )
        # y = (energy_1, energy_2, grad_1_proj)
        # b = np.linalg.lstsq(A, y, rcond=None)[0]
        # # np.poly1d expects the coefficients in decreasing power order
        # parabel = np.poly1d(b[::-1])
        # # Roots of the first derivative
        # roots = parabel.deriv().r
        # try:
            # real_minimum = (roots[roots.imag==0].real)[0]
        # except:
            # print("Only found an imaginary minimum for the fitted parabel. "
                  # "This should not happen :) Exiting.")
            # return

        # if abs(real_minimum) > 2*self.line_step_thresh:
            # print("Predicted minimum is above threshold! OHOHOHOHOH")
            # real_minimum = 0

        # self.mw_coords = mw_coords_1 + (real_minimum*D)
