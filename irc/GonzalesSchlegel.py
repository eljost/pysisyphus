#!/usr/bin/env python3

import logging
import warnings

import matplotlib.pyplot as plt
import numpy as np

#from pysisyphus.calculators.MullerBrownPot import MullerBrownPot
from pysisyphus.calculators.MullerBrownSympyPot import MullerBrownPot
from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.Geometry import Geometry

warnings.simplefilter('ignore', np.RankWarning)

# [1] http://pubs.acs.org/doi/pdf/10.1021/ja00295a002

def bfgs_update(H, grad_diffs, coord_diffs):
    y = grad_diffs[:,None]
    yT = grad_diffs[None,:]
    s = coord_diffs[:,None]
    sT = coord_diffs[None,:]
    #print("yTs", np.dot(yT, s))
    first_term = y.dot(yT) / yT.dot(s)
    second_term = H.dot(s).dot(sT).dot(H) / sT.dot(H).dot(s)
    dH = first_term - second_term
    return H + dH


def gs_step(geometry, max_step=0.15):
    summary = False
    hessian = geometry.hessian
    gradient0 = -geometry.forces
    gradient0_norm = np.linalg.norm(gradient0)
    # For the first BFGS update of the hessian we use differences
    # between the starting point of this micro-step and the initial
    # guess on the hypersphere.
    last_gradient = gradient0
    last_coords = geometry.coords

    # Determine pivot point x*_k+1.
    pivot_step = 0.5*max_step * gradient0/gradient0_norm

    # Take a step against the gradient to the pivot point.
    pivot_coords = geometry.coords - pivot_step

    # Make initial guess for x'_k+1.
    # For now just do a full step against the initial gradient.
    geometry.coords = pivot_coords - pivot_step
    # Initial displacement p' from the pivot point
    displacement = geometry.coords - pivot_coords
    print("initial displacement", np.linalg.norm(displacement))
    i = 0
    while True:
        #print(f"\tMicro: {i}")
        gradient = -geometry.forces
        gradient_diff = gradient - last_gradient
        coords_diff = geometry.coords - last_coords
        # After the first step move more or less along the surface of the
        # hypersphere. Now the BFGS updates are done between points on the
        # surface.
        last_gradient = gradient
        last_coords = geometry.coords
        # Update hessian and diagonalize it.
        #print(hessian)
        hessian = bfgs_update(hessian, gradient_diff, coords_diff)
        #print(hessian)
        eigvals, eigvecs = np.linalg.eig(hessian)
        #print("eigvecs", eigvecs)
        # Project gradient and displacement onto eigenvectors of
        # the hessian.
        gradients_proj = np.array([np.dot(gradient, v) for v in eigvecs])
        displacements_proj = np.array([np.dot(displacement, v) for v in eigvecs])
        # Find the lowest eigenvalue and use it to determine a first
        # guess for lambda.
        #lambda_ = eigvals[1]*.9
        lambda_ = np.sort(eigvals)[0] - 0.1
        prev_lambda = 0
        # Newton-Raphson to optimize lambda so f(lambda) = 0
        j = 0
        while abs(prev_lambda - lambda_) > 1e-10:
            # f(lambda)
            func = np.sum(
                    ((eigvals*displacements_proj-gradients_proj)
                     /(eigvals-lambda_))**2
            )
            func -= 0.25 * max_step**2
            # d(f(lambda))/(dlambda)
            deriv = np.sum(
                        2*(eigvals*displacements_proj-gradients_proj)**2
                        /(eigvals-lambda_)**3
            )
            prev_lambda = lambda_
            lambda_ = prev_lambda - func/deriv
            j += 1
            if j >= 100:
                logging.warning("Lambda search didn't converge!")
                break
        print(f"NR lambda search converged after {j} iterations: "
              f"λ={lambda_:5.2f}, f(λ)={func:8.5f}")
        # Calculate dx from optimized lambda
        dx = -np.dot(
                np.linalg.inv(hessian-lambda_*np.eye(hessian.shape[0])),
                gradient-lambda_*displacement
        )
        displacement += dx
        geometry.coords = pivot_coords + displacement
        #all_coords.append(geometry.coords)
        
        #break
        tangent = (gradient - (gradient.dot(displacement)
                              / np.linalg.norm(displacement)**2) * gradient
        )
        if (np.linalg.norm(dx) <= 1e-6):# and (np.linalg.norm(tangent) < 1e-6):
            print("WIN")
            summary = True
        if i >= 20:
            print("FAIL")
            summary = True
        if summary:
            print("new displ:", np.linalg.norm(displacement))
            print("constraint pTp - (1/2s)**2:",
                displacement[None,:].dot(displacement[:,None]) - 1/4*max_step**2
            )
            print("norm(dx)", np.linalg.norm(dx))
            print("norm(tangent)", np.linalg.norm(tangent))
            print()
            break

        i += 1

    return geometry.coords


def gonzales_schlegel(geometry, max_step=0.15):
    assert(max_step > 0), "max_step has to be > 0"

    all_coords = [geometry.coords, ]

    i = 0
    while i < 8:
        print(f"Macro: {i}")
        new_coords = gs_step(geometry, max_step)
        all_coords.append(new_coords)
        i += 1

    return np.array(all_coords)


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
    """
    print("all coordinates")
    for c in aic:
        print(c)
    """
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
