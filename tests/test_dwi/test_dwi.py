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
    true_grads = list()
    interpol_ens = list()
    interpol_grads = list()
    findiff = 1e-4
    fd_grads = list()
    for i in range(points):
        new_coords = coords[0] + i*step
        geom.coords = new_coords
        true_grad = geom.gradient
        true_grads.append(true_grad)
        true_energy = geom.energy
        true_ens.append(true_energy)
        interpol_energy, interpol_grad = dwi.interpolate(new_coords, gradient=True)
        interpol_ens.append(interpol_energy)
        interpol_grads.append(interpol_grad)

        fd_grad = list()
        for i in range(3):
            fd_step = np.zeros(3)
            fd_step[i] = findiff
            plus_val = dwi.interpolate(new_coords+fd_step)
            minus_val = dwi.interpolate(new_coords-fd_step)
            fd = (plus_val - minus_val) / (2*findiff)
            fd_grad.append(fd)
        fd_grads.append(fd_grad)

    true_ens = np.array(true_ens)
    interpol_ens = np.array(interpol_ens)
    true_grads = np.array(true_grads)
    interpol_grads = np.array(interpol_grads)
    fd_grads = np.array(fd_grads)
    np.testing.assert_allclose(interpol_grads, fd_grads)

    quiver_xs = np.arange(interpol_ens.size)
    quiver_true_ys = true_ens
    quiver_interpol_ys = interpol_ens
    fig, ax = plt.subplots()
    ax.plot(true_ens, label="True")
    ax.plot(interpol_ens, label="DWI")

    tg_U = true_grads[:,0]
    tg_V = true_grads[:,1]
    ti_U = interpol_grads[:,0]
    ti_V = interpol_grads[:,1]
    ax.quiver(quiver_xs, quiver_true_ys, tg_U, tg_V, color="b", label="grad(True)")
    ax.quiver(quiver_xs, quiver_interpol_ys, ti_U, ti_V, color="orange", label="grad(DWI)")
    ax.legend()

    plt.show()


if __name__ == "__main__":
    test_dwi()
