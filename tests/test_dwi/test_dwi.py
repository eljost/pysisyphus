#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.irc.DWI import DWI


def test_dwi():
    dwi = DWI()

    coords = np.array((
        (-0.222, 1.413, 0.),
        (-0.812, 1.242, 0.),
    ))
    geom = AnaPot.get_geom(coords[0])
    calc = geom.calculator

    c1 = geom.coords
    e1 = geom.energy
    g1 = geom.gradient
    h1 = geom.hessian
    dwi.update(c1, e1, g1, h1)

    geom.coords = coords[1]
    c2 = geom.coords
    e2 = geom.energy
    g2 = geom.gradient
    h2 = geom.hessian
    dwi.update(c2, e2, g2, h2)

    points = 10
    step = (coords[1] - coords[0]) / (points - 1)
    true_ens = list()
    interpol_ens = list()
    for i in range(points):
        new_coords = coords[0] + i*step
        geom.coords = new_coords
        true_energy = geom.energy
        true_ens.append(true_energy)
        interpol_energy = dwi.interpolate(new_coords)
        interpol_ens.append(interpol_energy)
    fig, ax = plt.subplots()
    ax.plot(true_ens, label="True")
    ax.plot(interpol_ens, label="DWI")
    ax.legend()

    plt.show()


if __name__ == "__main__":
    test_dwi()
