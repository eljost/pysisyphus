#!/usr/bin/env python3

import warnings

import matplotlib.pyplot as plt
import numpy as np

#from pysisyphus.calculators.MullerBrownPot import MullerBrownPot
from pysisyphus.calculators.MullerBrownSympyPot import MullerBrownPot
from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.Geometry import Geometry

warnings.simplefilter('ignore', np.RankWarning)

# [1] http://pubs.acs.org/doi/pdf/10.1021/ja00295a002

def bfgs_update(hessian, gradient_diff, coords_diff):
    return (hessian
                + np.dot(gradient_diff[:,None], gradient_diff[None,:])
                  / np.dot(gradient_diff[None,:], coords_diff[:,None])
                - np.dot(
                    np.dot(hessian, coords_diff[:,None]),
                    np.dot(coords_diff[None,:], hessian)
                )
                / np.dot(
                    np.dot(coords_diff[None,:], hessian),
                    coords_diff[:,None]
                )
                
    )


def gonzales_schlegel(geometry, max_step=0.15):
    assert(max_step > 0), "max_step has to be > 0"

    all_coords = [geometry.coords, ]
    all_gradients = list()

    # Determine pivot point x*_k+1.
    hessian = geometry.hessian
    gradient = -geometry.forces
    gradient_norm = np.linalg.norm(gradient)
    all_gradients.append(gradient)

    pivot_step = 0.5*max_step * gradient/gradient_norm
    pivot_coords = geometry.coords + pivot_step

    # Make initial guess for x'_k+1.
    # For now just do a full step along the initial gradient.
    geometry.coords = pivot_coords + pivot_step
    # p'
    displacement = geometry.coords - all_coords[-1]
    print("displacement")
    print(np.linalg.norm(displacement))
    i = 0
    while True:
        break
        print(f"Macro: {i}")
        cur_gradient = -geometry.forces
        gradient_diff = cur_gradient - all_gradients[-1]
        coords_diff = geometry.coords - all_coords[-1]
        # Update hessian and diagonalize it.
        hessian = bfgs_update(hessian, gradient_diff, coords_diff)
        eigvals, eigvecs = np.linalg.eig(hessian)
        # Project gradient and displacement onto eigenvectors of
        # the hessian.
        gradients_proj = np.array([np.dot(cur_gradient, v) for v in eigvecs])
        displacements_proj = np.array([np.dot(displacement, v) for v in eigvecs])
        # Find the lowest eigenvalue and use it to determine a first
        # guess for lambda.
        eigvals = np.sort(eigvals)
        lambda_ = eigvals[1]*.9
        #lambda_ = eigvals[0] - 0.1
        prev_lambda = 0
        # Newton-Raphson 
        j = 0
        while abs(prev_lambda - lambda_) > 1e-8:
            print(f"Micro {j}")
            # f(lambda)
            func = np.sum(
                    ((eigvals*displacements_proj-gradients_proj)
                     /(eigvals-lambda_))**2
            )
            func -= 0.25 * max_step**2
            # derivative of f(lambda)
            deriv = np.sum(
                        2*(eigvals*displacements_proj-gradients_proj)
                        /(eigvals-lambda_)**3
            )
            prev_lambda = lambda_
            lambda_ = prev_lambda - func/deriv
            print(f"lambda: {lambda_}")
            j += 1
        i += 1

    return np.array(all_coords)


    """
    i = 0
    while True:
        # Get gradient at Q0
        grad_0 = -geometry.forces
        energy_0 = geometry.energy

        if last_energy and (energy_0 > last_energy):
            print("Iteration {:04d}: Energy increased!".format(i))
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
        D_normed = D / np.linalg.norm(D)

        line_xs = [0, ]
        line_energies = [energy_1, ]

        line_step_size_thresh = 1.5*line_step_size
        # Try to find a useful point by projecting grad_1 on D.
        grad_1_normed = grad_1 / grad_1_norm
        step_D1 = np.dot(grad_1, D_normed) * D_normed * line_step_size
        step_D1_norm = np.linalg.norm(step_D1)
        if step_D1_norm < line_step_size_thresh:
            geometry.coords = coords_1 + step_D1
            step_D1_energy = geometry.energy
            line_xs.append(step_D1_norm)
            line_energies.append(step_D1_energy)
        # Otherwise just take a step along D.
        else:
            step_D2 = line_step_size * D_normed
            geometry.coords = coords_1 + step_D2
            step_D2_norm = np.linalg.norm(step_D2)
            line_xs.append(step_D2_norm)
            step_D2_energy = geometry.energy
            line_energies.append(step_D2_energy)

        # Calculate a 3rd point
        # Halve step size
        if line_energies[1] >= line_energies[0]:
            step_D3 = 0.5 * line_step_size * D_normed
        # Double step size
        else:
            step_D3 = 2 * line_step_size * D_normed

        step_D3_norm = np.linalg.norm(step_D3)
        geometry.coords = coords_1 + step_D3
        line_xs.append(step_D3_norm)
        step_D3_energy = geometry.energy
        line_energies.append(step_D3_energy)

        real_minimum = parabolic_fit(line_xs, line_energies)

        # ||step_norm|| = ||line_step_size * D_normed||
        #               = line_step_size * ||D_normed||
        #               = line_step_size * 1
        # So the resulting real_minimum basically corresponds
        # to a step size.
        irc_coords = coords_1 + (real_minimum*D_normed)
        geometry.coords = irc_coords

        last_energy = energy_0
        i += 1
        print()

    return np.array(all_irc_coords)
    """


def run():
    atoms = ("H", )

    # True TS
    calc, ts_coords = (MullerBrownPot(), np.array((-0.825, 0.624, 0.)))
    #calc, ts_coords = (MullerBrownPot(), np.array((-0.845041, 0.663752, 0.)))
    #calc, ts_coords = (MullerBrownPot(), np.array((-0.822, 0.624, 0.)))
    xlim = (-1.75, 1.25)
    ylim = (-0.5, 2.25)
    levels=(-150, -15, 40)

    #xlim = (-2, 2.5)
    #ylim = (0, 5)
    #levels = (-3, 4, 80)
    #calc, ts_coords = (AnaPot(), np.array((0.6906, 1.5491, 0.)))

    geometry = Geometry(atoms, ts_coords)
    geometry.set_calculator(calc)

    # Muller-Brown
    aic = gonzales_schlegel(geometry, max_step=0.15)
    # AnaPot
    #aic = gonzales_schlegel(geometry, max_step=0.0125)

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
