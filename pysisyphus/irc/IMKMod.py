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
            x_corr = (0, 0.5, 1)
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
        x = (0, self.line_step_size, factor*self.line_step_size)
        y = (energy_1, energy_2, energy_3)
        min_ = self.fit_parabola(x, y)
        step_corr = min_ * D
        mw_coords_4 = mw_coords_1 + step_corr
        self.mw_coords = mw_coords_4

        self.irc_coords.append(self.mw_coords)
        # self.irc_grads.append(self.gradient)
        self.irc_energies.append(self.energy)

        return
        # Project g1 along D
        grad_1_proj = grad_1.dot(D)
        """Fit parabel with 2 energies and 1 projected gradient.
        See [2] on how to do this. We want to determine the coefficients
        for the polynom f(x) = a0 + a1*x + a2*x**2. This can be done by
        solving A*b = y with b being the coefficient vector.
        The system matrix A is given by:
            A = (
                (1, x11, x11**2),
                (1, x12, x12**2),
                (0, 1, 2*x21)
            )
        energy_1 and and grad_1 are evaluated at the same point x11 = 0,
        so only x12 is != 0."""
        A = (
	    (1, 0, 0),
	    (1, self.line_step_size, self.line_step_size**2),
	    (0, 1, 0)
	)
        y = (energy_1, energy_2, grad_1_proj)
        b = np.linalg.lstsq(A, y, rcond=None)[0]
        # np.poly1d expects the coefficients in decreasing power order
        parabel = np.poly1d(b[::-1])
        # Roots of the first derivative
        roots = parabel.deriv().r
        try:
            real_minimum = (roots[roots.imag==0].real)[0]
        except:
            print("Only found an imaginary minimum for the fitted parabel. "
                  "This should not happen :) Exiting.")
            return

        if abs(real_minimum) > 2*self.line_step_thresh:
            print("Predicted minimum is above threshold! OHOHOHOHOH")
            real_minimum = 0

        self.mw_coords = mw_coords_1 + (real_minimum*D)


def run():
    from pysisyphus.calculators.MullerBrownPot import MullerBrownPot
    from pysisyphus.calculators.AnaPot import AnaPot
    atoms = ("H", )

    # True TS
    # calc, ts_coords = (MullerBrownPot(), np.array((-0.825, 0.624, 0.)))
    #calc, ts_coords = (MullerBrownPot(), np.array((-0.845041, 0.663752, 0.)))
    # xlim = (-1.75, 1.25)
    # ylim = (-0.5, 2.25)
    # levels=(-150, -15, 40)

    xlim = (-2, 2.5)
    ylim = (0, 5)
    levels = (-3, 4, 80)
    calc, ts_coords = (AnaPot(), np.array((0.61173, 1.49297, 0.)))

    geometry = Geometry(atoms, ts_coords)
    geometry.set_calculator(calc)

    # Muller-Brown
    #aic = imk_mod(geometry, desired_step=0.15, line_step_size=0.05)
    # aic = imk_mod(geometry, desired_step=0.3)
    # AnaPot
    aic = IMKMod(geometry, step_length=0.2, max_steps=50)#, backward=False)
    from pysisyphus.irc.Euler import Euler
    # aic = Euler(geometry, step_length=0.2, max_steps=75)
    aic.run()
    
    fig, ax = plt.subplots(figsize=(8,8))

    # Calculate the potential
    x = np.linspace(*xlim, 100)
    y = np.linspace(*ylim, 100)
    X, Y = np.meshgrid(x, y)
    Z = np.full_like(X, 0)
    fake_atoms = ("H", )
    pot_coords = np.stack((X, Y, Z))
    pot = calc.get_energy(fake_atoms, pot_coords)["energy"]
    levels = np.linspace(*levels)
    contours = ax.contour(X, Y, pot, levels)

    ax.scatter(*ts_coords[:2], s=75, c="k")
    ax.scatter(-1.05274, 1.02776, s=75, c="k")
    ax.scatter(1.94101, 3.85427, s=75, c="k")

    def label_steps(ax, xs, ys):
        for i, (x, y) in enumerate(zip(xs, ys)):
            ax.annotate(f"{i}", (x, y))

    try:
        fw_coords = np.array(aic.forward_coords)
        ax.plot(fw_coords[:, 0], fw_coords[:, 1], "ro", ls="-", label="forward")
        label_steps(ax, fw_coords[:,0][::-1], fw_coords[:,1][::-1])
    except AttributeError:
        pass
    try:
        bw_coords = np.array(aic.backward_coords)
        ax.plot(bw_coords[:, 0], bw_coords[:, 1], "bo", ls="-", label="backward")
        label_steps(ax, bw_coords[:,0], bw_coords[:,1])
    except AttributeError:
        pass
    ax.legend()
    plt.show()


if __name__ == "__main__":
    run()
