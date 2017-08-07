#!/usr/bin/env python3

import warnings

import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.MullerBrownPot import MullerBrownPot
from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.Geometry import Geometry

warnings.simplefilter('ignore', np.RankWarning)

# [1] http://pubs.acs.org/doi/pdf/10.1021/ja00295a002


def parabolic_fit(xs, ys):
    fit = np.polyfit(xs, ys, deg=2)
    fit = np.poly1d(fit)
    minima = fit.deriv().r
    real_minima = minima[minima.imag==0].real

    return real_minima


def imk_mod(geometry, desired_step=0.15, line_step_size=0.025, fit_points=3):
    all_irc_coords = list()

    last_energy = None
    i = 0
    while True:
        # Get gradient at Q0
        grad_0 = -geometry.forces
        energy_0 = geometry.energy

        if last_energy and (energy_0 > last_energy):
            print("Iteration {:04d}: Energy incresed!".format(i))
            break

        all_irc_coords.append(geometry.coords)

        grad_0_norm = np.linalg.norm(grad_0)
        step_size = np.linalg.norm(desired_step) / grad_0_norm
        step = -step_size * grad_0

        # Step to intermediate point Q1
        coords_1 = geometry.coords + step
        geometry.coords = coords_1

        # Get gradient at Q1
        grad_1 = -geometry.forces
        energy_1 = geometry.energy
        grad_1_norm = np.linalg.norm(grad_1)

        # Determine bisector
        D = grad_0/grad_0_norm - grad_1/grad_1_norm

        line_xs = [0, ]
        line_energies = [energy_1, ]

        line_step_size_thresh = 2*line_step_size
        # Try to find a useful point by projecting grad_1 on D.
        step_D1 = np.dot(grad_1, D) * D * line_step_size
        step_D1_norm = np.linalg.norm(step_D1)
        if step_D1_norm < line_step_size_thresh:
            geometry.coords = coords_1 + step_D1
            step_D1_energy = geometry.energy
            line_xs.append(step_D1_norm)
            line_energies.append(step_D1_energy)
        # Otherwise just take a step along D.
        else:
            step_D2 = line_step_size * D
            geometry.coords = coords_1 + step_D2
            step_D2_norm = np.linalg.norm(step_D2)
            line_xs.append(step_D2_norm)
            step_D2_energy = geometry.energy
            line_energies.append(step_D2_energy)

        real_minimum = parabolic_fit(line_xs, line_energies)
        print("real_minimum", real_minimum)
        #if abs(real_minimum) > abs(line_step_size_thresh):

        """
        line_steps = [j * D for j in range(3)]
        line_energies = list()
        for line_step in line_steps:
            line_search_coords = coords_1 + line_step
            geometry.coords = line_search_coords
            energy = geometry.energy
            line_energies.append(energy)
        line_xs = np.linalg.norm(line_steps, axis=-1)
        line_energies = np.array(line_energies)

        real_minimum = parabolic_fit(line_xs, line_energies)
        """

        irc_coords = coords_1 + (real_minimum*line_step_size*D)
        geometry.coords = irc_coords

        last_energy = energy_0
        i += 1

    return np.array(all_irc_coords)


def run():
    #ts_coords = np.array((-0.822, 0.624, 0.))
    atoms = ("H", )
    calc, ts_coords = (MullerBrownPot(), np.array((-0.845041, 0.663752, 0.)))
    xlim = (-1.75, 1.25)
    ylim = (-0.5, 2.25)
    levels=(-150, -15, 40)

    #xlim = (-2, 2.5)
    #ylim = (0, 5)
    #levels = (-3, 4, 80)
    #calc, ts_coords = (AnaPot(), np.array((0.6906, 1.5491, 0.)))
    geometry = Geometry(atoms, ts_coords)
    geometry.set_calculator(calc)

    aic = imk_mod(geometry)

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

    ax.plot(aic[:, 0], aic[:, 1], "ro", ls="-")
    plt.show()


if __name__ == "__main__":
    run()
