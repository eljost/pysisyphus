#!/usr/bin/env python3

import logging
import warnings

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import newton

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


def bracket(func, a, b):
    max_cycles = 50
    factor = 1.6
    fa = func(a)
    fb = func(b)
    for i in range(max_cycles):
        if fa*fb < 0.0:
            break
        if (abs(fa) < abs(fb)):
            a += factor * (a-b)
            fa = func(a)
        else:
            b += factor * (b-a)
            fb = func(b)
    if i > max_cycles:
        a = None
        b = None
    return a, b


def get_lower_bracket(func, lower, upper):
    max_cycles = 50
    factor = 1.5
    f_lower = func(lower)
    f_upper = func(upper)
    while (f_lower * f_upper) > 0.0:
        assert(lower-upper) < 0
        lower += factor * (lower-upper)
        f_lower = func(lower)

    return lower

    """
    for i in range(max_cycles):
        if f_lower*f_upper< 0.0:
            break
        if (abs(fa) < abs(fb)):
            a += factor * (a-b)
            fa = func(a)
        else:
            b += factor * (b-a)
            fb = func(b)
    if i > max_cycles:
        a = None
        b = None
    return a, b
    """



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


def gs_step(geometry, max_step):
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
    last_lambda = None
    eye = np.eye(displacement.size)
    while True:
        print(f"Micro {i}")
        gradient = -geometry.forces
        gradient_diff = gradient - last_gradient
        coords_diff = geometry.coords - last_coords
        # After the first step move more or less along the surface of the
        # hypersphere. Now the BFGS updates are done between points on the
        # surface.
        last_gradient = gradient
        last_coords = geometry.coords
        #hessian = bfgs_update(hessian, gradient_diff, coords_diff)
        hessian = geometry.hessian
        eigvals, eigvecs = np.linalg.eig(hessian)
        hinv = np.linalg.inv(hessian)

        def ffunc(lambda_):
            hmlinv = np.linalg.inv(hessian - eye*lambda_)
            glp = gradient - displacement*lambda_
            #print("hmlinv", hmlinv.dot(hessian - eye*lambda_))
            tmp = displacement - hmlinv.dot(glp)
            return tmp.dot(tmp) - 0.25*(max_step**2)

        smallest_eigval = np.sort(eigvals)[0]
        lambda_ = np.sort(eigvals)[0]
        lambda_ *= 1.5 if (lambda_ < 0) else 0.5
        lambda_guess = lambda_
        lower_bracket = get_lower_bracket(ffunc, lambda_, smallest_eigval)
        newton_lambda = newton(ffunc, lambda_)
        lambda_ = newton_lambda
        print("lower_bracket", lower_bracket)
        print("newton_lambda", newton_lambda)
        """
        if abs(lambda_guess - lambda_) < 50:
            print("small!")
            lambda_ -= 1500
        """
        print("lambda_guess", lambda_, "eigvals", eigvals)
        #print(hessian)
        #print(hinv)
        #print(hessian.dot(hinv))

        """
        lambdas = np.linspace(-10000, 10000, 100000)
        plt.plot(lambdas, [ffunc(l) for l in lambdas])
        plt.show()
        import sys; sys.exit()
        """


        # Project gradient and displacement onto eigenvectors of
        # the hessian.
        gradients_proj = np.array([np.dot(gradient, v) for v in eigvecs])
        displacements_proj = np.array([np.dot(displacement, v) for v in eigvecs])

        def myffunc(lambda_):
            f = np.sum(
                    ((eigvals*displacements_proj-gradients_proj)
                     /(eigvals-lambda_))**2
            )
            return f- 0.25 * max_step**2
        my_newton_lambda = newton(myffunc, lambda_guess)
        print("my_newton_lambda", my_newton_lambda)

        """
        # Find the lowest eigenvalue and use it to determine a first
        # guess for lambda.
        lambda_ = np.sort(eigvals)[0] - 0.1
        bracket_guess = (2*lambda_, lambda_)
        low_brack = min(bracket_guess)
        up_brack = max(bracket_guess)
        low_brack, up_brack = bracket(ffunc, low_brack, up_brack)
        #lambdas = np.linspace(low_brack, up_brack, 100)
        print("low_brack", low_brack, "up_brack", up_brack)
        if not last_lambda:
            last_lambda = lambda_
        """

        """
        prev_lambda = 0
        # Newton-Raphson to optimize lambda so f(lambda) = 0
        j = 0
        while abs(prev_lambda - lambda_) > 1e-12:
            prev_lambda = lambda_
            # f(lambda)
            func = np.sum(
                    ((eigvals*displacements_proj-gradients_proj)
                     /(eigvals-lambda_))**2
            )
            func -= 0.25 * max_step**2
            # d(f(lambda))/(dlambda)
            deriv = 2*np.sum(
                        (eigvals*displacements_proj-gradients_proj)**2
                        /(eigvals-lambda_)**3
            )
            lambda_ = prev_lambda - func/deriv
            j += 1
            if j >= 100:
                logging.warning("Lambda search didn't converge!")
                break
        
        print(f"NR lambda search converged after {j} iterations: "
              f"λ={lambda_:5.5f}, f(λ)={func:8.10f}")
        """

        """
        lambdas = np.linspace(-3000, 3000, 10000)
        plt.plot(lambdas, [ffunc(l) for l in lambdas])
        ev_fmt = ", ".join(["{:.1f}".format(ev) for ev in eigvals])
        plt.title("ev {}, λ {:.1f}".format(ev_fmt, lambda_guess))
        plt.ylim(-0.05, 0.05)
        plt.axhline(0)
        plt.show()
        """

        """
        print("last_lambda", last_lambda)
        print("lambda_", lambda_)
        print("l / ll", lambda_ / last_lambda)
        if (lambda_ / last_lambda) > 5:
            lambda_ = last_lambda
            print("HMM!")
        else:
            last_lambda = lambda_
        """
        # Calculate dx from optimized lambda
        dx = -np.dot(
                np.linalg.inv(hessian-lambda_*np.eye(hessian.shape[0])),
                gradient-lambda_*displacement
        )
        print("norm(dx)", np.linalg.norm(dx))
        displacement += dx
        print("norm(displacement)", np.linalg.norm(displacement))
        geometry.coords = pivot_coords + displacement
        micro_coords.append(geometry.coords)
        
        if (np.linalg.norm(dx) <= 1e-5):
            print("WIN")
            summary = True
        elif i >= 5:
            print("FAIL")
            summary = True
            #import sys; sys.exit()

        if summary:
            print("new displ:", np.linalg.norm(displacement))
            print("constraint pTp - (1/2s)**2:",
                displacement.dot(displacement) - 1/4*max_step**2
            )
            print("norm(dx)", np.linalg.norm(dx))
            print()
            displ_norm = np.linalg.norm(displacement)
            tangent = gradient - gradient.dot(displacement)/displ_norm * gradient
            """
            tangent = (gradient - (gradient.dot(displacement)
                                  / (np.linalg.norm(displacement)**2)) * gradient
            )
            """
            print("tangent:", tangent, "norm(tangent)", np.linalg.norm(tangent))
            break

        i += 1
        print()

    return geometry.coords, pivot_coords, init_displ, micro_coords


def gonzales_schlegel(geometry, max_step=0.15):
    assert(max_step > 0), "max_step has to be > 0"

    all_coords = [geometry.coords, ]
    pivot_coords = list()
    init_guess_coords = list()
    micro_coords_list = list()


    i = 0
    while i < 10:
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
    calc, ts_coords = (MullerBrownSympyPot2D(), np.array((-0.845041, 0.663752)))
    #calc, ts_coords = (MullerBrownSympyPot2D(), np.array((-0.98137551,  0.91393777)))
    #calc, ts_coords = (MullerBrownPot(), np.array((-0.822, 0.624, 0.)))
    #xlim = (-1.75, 1.25)
    #ylim = (-0.5, 2.25)
    xlim = (-1.25, -.25)
    ylim = (0.5, 1.5)
    levels=(-150, -15, 40)

    #xlim = (-2, 2.5)
    #ylim = (0, 5)
    #levels = (-3, 4, 80)
    #calc, ts_coords = (AnaPot(), np.array((0.6906, 1.5491, 0.)))
    #calc, ts_coords = (AnaPot2D(), np.array((0.6906, 1.5491)))
    #calc, ts_coords = (AnaPot2D(), np.array((0.830, 1.67)))

    geometry = Geometry(atoms, ts_coords)
    geometry.set_calculator(calc)

    # Muller-Brown
    aic, pc, igc, mcl = gonzales_schlegel(geometry, max_step=0.1)
    print("all coordinates")
    for i, c in enumerate(aic):
        print(i, c)
    # AnaPot2D
    #aic, pc, igc, mcl = gonzales_schlegel(geometry, max_step=0.2)

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
        for i, m in enumerate(mc):
            ax.text(*m, str(i))
    ax.plot(aic[:, 0], aic[:, 1], "ro", ls="-")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    run()
