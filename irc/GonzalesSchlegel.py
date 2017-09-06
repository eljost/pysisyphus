#!/usr/bin/env python3

import logging
import warnings

import matplotlib.pyplot as plt
import numpy as np

#from pysisyphus.calculators.MullerBrownSympyPot import MullerBrownPot
from pysisyphus.calculators.MullerBrownSympyPot2D import MullerBrownSympyPot2D
from pysisyphus.calculators.AnaPot2D import AnaPot2D
from pysisyphus.Geometry import Geometry


def bfgs_update(H, grad_diffs, coord_diffs):
    y = grad_diffs[:,None]
    yT = grad_diffs[None,:]
    s = coord_diffs[:,None]
    sT = coord_diffs[None,:]
    first_term = y.dot(yT) / yT.dot(s)
    second_term = H.dot(s).dot(sT).dot(H) / sT.dot(H).dot(s)
    dH = first_term - second_term
    return H + dH

def bfgs_update_ase(H, grad_diffs, coord_diffs):
    df = -grad_diffs
    dr = coord_diffs
    a = np.dot(dr, df)
    dg = np.dot(H, dr)
    b = np.dot(dr, dg)
    dH = -(np.outer(df, df) / a + np.outer(dg, dg) / b)
    print("dH ase")
    print(dH)
    return H + dH


def newton_raphson(f, df, a, b, tol=1e-8):
    fa = f(a)
    if fa == 0.0:
        return a
    fb = f(b)
    if fb == 0.0:
        return b
    if fa*fb > 0.0:
        raise Exception("Root is not bracketed!")
    x = 0.5 * (a+b)
    for i in range(30):
        fx = f(x)
        if abs(fx) < tol: return x
        # Tighten the brackets on the root
        if fa*fx < 0.0:
            b = x
        else:
            a = x
        # Try a Newton-Raphson step
        dfx = df(x)
        # If division by zero push out of bound
        try:
            dx = -fx/dfx
        except ZeroDivisionError:
            dx = b - a
        x = x + dx
        # If the result is outside the brackets, use bisection
        if (b-x)*(x-a) < 0.0:
            dx = 0.5*(b-a)
            x = a + dx
        # Check for convergence
        if abs(dx) < tol*max(abs(b),1.0):
            return x



def gs_step(geometry, max_step=0.15):
    summary = False
    micro_coords = list()
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
    init_displ = displacement
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
        #hessian = bfgs_update(hessian, gradient_diff, coords_diff)
        #print("updated hessian")
        #print(hessian)
        anal_hessian = geometry.hessian
        #print("analytical hessian")
        #print(anal_hessian)
        #print("diff")
        #print(anal_hessian-hessian)
        #print()
        hessian = anal_hessian
        eigvals, eigvecs = np.linalg.eig(hessian)
        #print("eigvecs", eigvecs)

        # Project gradient and displacement onto eigenvectors of
        # the hessian.
        gradients_proj = np.array([np.dot(gradient, v) for v in eigvecs])
        displacements_proj = np.array([np.dot(displacement, v) for v in eigvecs])
        # Find the lowest eigenvalue and use it to determine a first
        # guess for lambda.
        lambda_ = np.sort(eigvals)[0] - 0.1

        def func(lambda_):
            f = np.sum(
                    ((eigvals*displacements_proj-gradients_proj)
                     /(eigvals-lambda_))**2
            )
            f -= 0.25 * max_step**2
            return f

        def deriv(lambda_):
            return np.sum(
                        2*(eigvals*displacements_proj-gradients_proj)**2
                        /(eigvals-lambda_)**3
            )
        #a, b = (-10, 10)
        #lambda_ = newton_raphson(func, deriv, a, b, tol=1e-8)
        """
        print("diff", abs(lambda_ + 0.1))
        if abs(lambda_ + 0.1) <= 1e-4:
            print("EY")
            lambda_ = eigvals[1]-0.1
        """
        print("lambda guess", lambda_)
        prev_lambda = 0
        j = 0
        lambda_lb = -10
        lambda_ub = 10
        # Newton-Raphson to optimize lambda so f(lambda) = 0
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
            #print(func)
            j += 1
            if j >= 200:
                logging.warning("Lambda search didn't converge!")
                break
        """
        while True:
            func = 0
            df = 0
            for k in range(eigvals.size):
                bpgb = (eigvals[k]*displacements_proj[k]-gradients_proj[k])/(eigvals[k]-lambda_)
                bl = eigvals[k]- lambda_
                func += bpgb * bpgb
                df += 2.0 * bpgb*bpgb / bl
            func -= 0.25*max_step*max_step
            prev_lambda = lambda_
            lambda_ = lambda_ - func/df
            #print(lambda_)
            if abs(prev_lambda - lambda_) < 1e-8:
                break
        """
        print(f"NR lambda search converged after {j} iterations: "
              f"λ={lambda_:5.5f}, f(λ)={func:8.10f}")
        # Calculate dx from optimized lambda
        dx = -np.dot(
                np.linalg.inv(hessian-lambda_*np.eye(hessian.shape[0])),
                gradient-lambda_*displacement
        )
        print("dx", dx)
        displacement += dx
        geometry.coords = pivot_coords + displacement
        micro_coords.append(geometry.coords)
        #all_coords.append(geometry.coords)
        
        #break
        tangent = (gradient - (gradient.dot(displacement)
                              / np.linalg.norm(displacement)**2) * gradient
        )
        if (np.linalg.norm(dx) <= 1e-5):# and (np.linalg.norm(tangent) < 1e-6):
            print("WIN")
            summary = True
        elif i >= 20:
            print("FAIL")
            summary = True
            #import sys; sys.exit()

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

    return geometry.coords, pivot_coords, init_displ, micro_coords


def gonzales_schlegel(geometry, max_step=0.15):
    assert(max_step > 0), "max_step has to be > 0"

    all_coords = [geometry.coords, ]
    pivot_coords = list()
    init_guess_coords = list()
    micro_coords_list = list()


    i = 0
    while i < 1:
        print(f"Macro: {i}")
        new_coords, pv, init_displ, micro_coords = gs_step(geometry, max_step)
        all_coords.append(new_coords)
        pivot_coords.append(pv)
        init_guess_coords.append(pv+init_displ)
        micro_coords_list.append(micro_coords)
        i += 1

    return (np.array(all_coords), np.array(pivot_coords),
            np.array(init_guess_coords), np.array(micro_coords_list))


def run():
    atoms = ("H", )

    # True TS
    #calc, ts_coords = (MullerBrownPot(), np.array((-0.825, 0.624, 0.)))
    #calc, ts_coords = (MullerBrownPot(), np.array((-0.845041, 0.663752, 0.)))
    #calc, ts_coords = (MullerBrownSympyPot2D(), np.array((-0.845041, 0.663752)))
    #calc, ts_coords = (MullerBrownPot(), np.array((-0.822, 0.624, 0.)))
    #xlim = (-1.75, 1.25)
    #ylim = (-0.5, 2.25)
    #xlim = (-1.25, -.25)
    #ylim = (0.5, 1.5)
    #levels=(-150, -15, 40)

    xlim = (-2, 2.5)
    ylim = (0, 5)
    levels = (-3, 4, 80)
    #calc, ts_coords = (AnaPot(), np.array((0.6906, 1.5491, 0.)))
    calc, ts_coords = (AnaPot2D(), np.array((0.6906, 1.5491)))

    geometry = Geometry(atoms, ts_coords)
    geometry.set_calculator(calc)

    # Muller-Brown
    #aic, pc, igc, mcl = gonzales_schlegel(geometry, max_step=0.15)
    """
    print("all coordinates")
    for c in aic:
        print(c)
    """
    # AnaPot2D
    aic, pc, igc, mcl = gonzales_schlegel(geometry, max_step=0.2)

    fig, ax = plt.subplots(figsize=(8,8))

    """
    # Calculate the potential
    # 3D
    x = np.linspace(*xlim, 100)
    y = np.linspace(*ylim, 100)
    X, Y = np.meshgrid(x, y)
    Z = np.full_like(X, 0)
    fake_atoms = ("H", )
    pot_coords = np.stack((X, Y, Z))
    """
    # 2D
    x = np.linspace(*xlim, 100)
    y = np.linspace(*ylim, 100)
    X, Y = np.meshgrid(x, y)
    fake_atoms = ("H", )
    pot_coords = np.stack((X, Y))
    pot = calc.get_energy(fake_atoms, pot_coords)["energy"]
    levels = np.linspace(*levels)
    contours = ax.contour(X, Y, pot, levels)

    ax.plot(pc[:, 0], pc[:, 1], "bo", ls="-", label="pivot")
    ax.plot(igc[:, 0], igc[:, 1], "go", ls="-", label="initial guess")
    for mc in mcl:
        mc = np.array(mc)
        ax.plot(mc[:, 0], mc[:, 1], "yo", ls="-")
    ax.plot(aic[:, 0], aic[:, 1], "ro", ls="-")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    run()
