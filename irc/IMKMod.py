#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.MullerBrownPot import MullerBrownPot
from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.Geometry import Geometry

# http://pubs.acs.org/doi/pdf/10.1021/ja00295a002

def imk_irc(geometry, desired_step=0.5, line_step_size=0.025):
    all_irc_coords = [geometry.coords.copy(), ]

    for i in range(15):
        print(f"ITER {i}")
        # Get gradient at Q0
        grad_0 = -geometry.forces
        grad_0_norm = np.linalg.norm(grad_0)
        step_size = np.linalg.norm(desired_step) / grad_0_norm
        step = -step_size * grad_0
        print("step", step)
        coords_1 = geometry.coords + step
        #print(step)
        # Step to intermediate point
        geometry.coords = coords_1
        # Get gradient at Q1
        grad_1 = -geometry.forces
        print("grad_1", grad_1)

        # Determine bisector
        print("g0_norm", grad_0_norm)
        grad_1_norm = np.linalg.norm(grad_1)
        print("g1_norm", grad_1_norm)
        D = grad_0/grad_0_norm - grad_1/grad_1_norm
        print("D", D)

        js = np.arange(7)
        line_energies = list()
        # Store Q1 to be used in the linesearch
        coords_backup = geometry.coords.copy()
        for j in js:
            line_step = j*line_step_size*D
            line_search_coords = coords_backup + line_step
            geometry.coords = line_search_coords
            energy = geometry.energy
            print("line_step", line_step)
            print("j", j, "energy", energy)
            line_energies.append(energy)
        line_energies = np.array(line_energies)
        fit = np.polyfit(js, line_energies, deg=2)
        fit = np.poly1d(fit)
        print("fit")
        print(fit)
        crit = fit.deriv().r
        print("crit", crit)
        r_crit = crit[crit.imag==0].real
        print("r_crit", r_crit)
        test = fit.deriv(2)(r_crit) 
        print("test", test)
        irc_coords = coords_backup + (r_crit[0]*line_step_size*D)
        all_irc_coords.append(irc_coords)
        geometry.coords = irc_coords

        print("New_step", step_size)

        print()

    return np.array(all_irc_coords)


def run():
    #ts_coords = np.array((-0.822, 0.624, 0.))
    atoms = ("H", )
    #calc, ts_coords = (MullerBrownPot(), np.array((-0.845041, 0.663752, 0.)))
    #xlim = (-1.75, 1.25)
    #ylim = (-0.5, 2.25)
    #levels=(-150, -15, 40)

    xlim = (-2, 2.5)
    ylim = (0, 5)
    levels = (-3, 4, 80)
    calc, ts_coords = (AnaPot(), np.array((0.6906, 1.5491, 0.)))
    geometry = Geometry(atoms, ts_coords)
    geometry.set_calculator(calc)

    aic = imk_irc(geometry)

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
