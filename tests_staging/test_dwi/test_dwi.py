#!/usr/bin/env python3

from pprint import pprint

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


def test_euler():
    dwi = DWI()

    coords = np.array((
        (-0.222, 1.413, 0.),
        (-0.812, 1.242, 0.),
    ))
    # Half step
    # diff = coords[1] - coords[0]
    # coords[1] = coords[0] + 0.5*diff
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

    # Euler integration
    norm = np.linalg.norm(coords[1] - coords[0])
    print(f"             norm={norm:.8f}")
    all_coords = list()
    richardson = dict()
    errors = list()
    for k in range(10):
        points = 10*(2**k) + 1
        corr_step_length  = norm / (points - 1)
        # print("corr_step_length", corr_step_length)
        cur_coords = coords[0].copy()
        k_coords = list()
        length = 0
        while True:
            k_coords.append(cur_coords.copy())
            if length >= norm:
                # print(f"Converged! length={length:.8f}, length-step={length-corr_step_length:.8f}")
                print(f"Converged! length={length:.8f}")
                break
            energy, gradient = dwi.interpolate(cur_coords, gradient=True)
            cur_coords += corr_step_length * -gradient/np.linalg.norm(gradient)
            length += corr_step_length
            # Check for oscillation
            try:
                prev_coords = k_coords[-2]
                osc_norm = np.linalg.norm(cur_coords - prev_coords)
                if osc_norm <= corr_step_length:
                    print("Detected oscillation. Breaking!")
                    # TODO: handle this by restarting everyhting with a smaller stepsize.
                    # Check 10.1039/c7cp03722h SI
                    assert False, "This case is not yet handled"
                    break
            except IndexError:
                pass
        richardson[(k, 0)] = cur_coords

        # Refine using Richardson extrapolation
        # Set additional values using Richard extrapolation
        for j in range(1, k+1):
            print(f"k={k},j={j}")
            richardson[(k, j)] = ((2**j) * richardson[(k, j-1)] - richardson[(k-1, j-1)]) \
                                 / (2**j-1)
        if k > 0:
            # RMS of coordinates
            error = np.sqrt(np.mean((richardson[(k, k)] - richardson[(k-1, k-1)])**2))
            print(f"\terror={error:.8e}")
            errors.append(error)
            if error <= 1e-6:
                break
        all_coords.append(np.array(k_coords))

    print("Richardson table")
    pprint(richardson)
    print()
    print("Errors")
    erros = np.array(errors)
    print(errors)

    calc = geom.calculator
    calc.plot()
    ax = calc.ax
    ax.plot(*coords.T[:2], "r-")

    for i, ac in enumerate(all_coords, 1):
        ax.plot(*ac.T[:2], label=i)
    ax.legend()
    plt.show()


if __name__ == "__main__":
    # test_dwi()
    test_euler()
