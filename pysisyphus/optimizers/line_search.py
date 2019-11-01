#!/usr/bin/env python3

# Based on
# [1] Nocedal - Numerical Optimization, chapter 3. Line Search Methods

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import line_search as sp_line_search


def interpol_alpha_quad(f_0, df_0, f_alpha_0, alpha_0):
    return -df_0*alpha_0**2 / 2 / (f_alpha_0 - f_0 - df_0*alpha_0)


def interpol_alpha_cubic(f_0, df_0, f_alpha_0, f_alpha_1, alpha_0, alpha_1):
    quot = 1 / (alpha_0**2 * alpha_1**2 * (alpha_1 - alpha_0))
    A = np.array(((alpha_0**2, -alpha_1**2),
                  (-alpha_0**3, alpha_1**3)))
    B = np.array(((f_alpha_1 - f_0 - df_0*alpha_1),
                  (f_alpha_0 - f_0 - df_0*alpha_0)))
    a, b = quot * A @ B
    alpha_cubic = (-b + (b**2 - 3*a*df_0)**(0.5)) / (3*a)
    return alpha_cubic


def line_search(f, df, x, p, f_0=None, f_0_prev=None, df_0=None,
                alpha=None, alpha_min=0.01, alpha_max=100.,
                c1=1e-4, c2=0.9, fac=2, strong=True):

    alpha_fs = {}
    def phi_alpha(alpha):
        try:
            f_alpha = alpha_fs[alpha]
        except KeyError:
            f_alpha = f(x + alpha*p)
            alpha_fs[alpha] = f_alpha
        return f_alpha

    alpha_dfs = {}
    def dphi_alpha(alpha):
        try:
            df_alpha = alpha_dfs[alpha]
        except KeyError:
            df_alpha = df(x + alpha*p)
            alpha_dfs[alpha] = df_alpha
        return df_alpha@p

    def results_for(alpha):
        return alpha, alpha_fs[alpha], alpha_dfs[alpha]

    # We shouldn't have to recompute phi_0 and dphi_0, as they were probably already
    # calculated when determining the descent direction p.
    if f_0 is None:
        phi_0 = phi_alpha(0)
    else:
        phi_0 = f_0
        alpha_fs[0] = f_0
    if df_0 is None:
        dphi_0 = dphi_alpha(0)
    else:
        dphi_0 = df_0 @ p
        alpha_dfs[0] = df_0

    def sufficiently_decreased(phi_alpha, alpha):
        return phi_alpha <= (phi_0 + c1 * alpha * dphi_0)

    def curvature_condition(dphi_alpha):
        return dphi_alpha >= c2*dphi_0

    def curvature_condition_strong(dphi_alpha):
        return abs(dphi_alpha) <= -c2*dphi_0

    curv_cond = curvature_condition_strong if strong else curvature_condition

    def zoom(alpha_lo, alpha_hi, phi_0, dphi_0, phi_lo,
             phi_alpha_=None, alpha_0_=None, max_cycles=10):

        alphas = list()
        phi_alphas = list()
        if phi_alpha_:
            phi_alphas = [phi_alpha_, ]
        if alpha_0_:
            alphas = [alpha_0_, ]

        for j in range(max_cycles):
            # Interpoaltion of alpha between alpha_lo, alpha_hi
            #
            # Try cubic interpolation if at least two additional alphas and
            # corresponding phi_alpha values are available beside alpha = 0.
            if len(phi_alphas) > 1:
                alpha_prev = alphas[-1]
                phi_alpha_prev = phi_alphas[-1]
                alpha_j = interpol_alpha_cubic(phi_0, dphi_0,
                                               phi_alpha_, phi_alpha_prev,
                                               alpha_0_, alpha_prev
                )
            # Try quadratic interpolation if at one additional alpha and
            # corresponding phi_alpha value is available beside alpha = 0.
            elif len(phi_alphas) == 1:
                alpha_j = interpol_alpha_quad(phi_0, dphi_0, phi_alpha_, alpha_0_)
            # Fallback to simple bisection
            else:
                alpha_j = (alpha_lo + alpha_hi) / 2

            phi_j = phi_alpha(alpha_j)
            # Store the values so they can be reused for cubic interpolation
            alphas.append(alpha_j)
            phi_alphas.append(phi_j)

            # True if alpha is still too big or if the function value
            # increased compared to the previous cycle.
            if (not sufficiently_decreased(phi_j, alpha_j)
                or phi_j > phi_lo):
                # Shrink interval to (alpha_lo, alpha_j)
                alpha_hi = alpha_j
                continue

            dphi_j = dphi_alpha(alpha_j)
            if curv_cond(dphi_j):
                print(f"\tzoom converged after {j+1} cycles.")
                return alpha_j

            if (dphi_j * (alpha_hi - alpha_lo)) >= 0:
                alpha_hi = alpha_lo
            # Shrink interval to (alpha_j, alpha_hi)
            alpha_lo = alpha_j
        raise Exception("zoom() didn't converge in {j+1} cycles!")

    
    alpha_prev = 0
    phi_prev = phi_0
    if alpha is not None:
        alpha_i = alpha
    # This does not seem to help
    # elif f_0_prev is not None:
        # alpha_i = min(1.01*2*(f_0 - f_0_prev) / dphi_0, 1.)
        # print("ai", alpha_i)
        # alpha_i = 1. if alpha_i < 0. else alpha_i
    else:
        alpha_i = 1.0
    for i in range(50):
        phi_i = phi_alpha(alpha_i)
        if (not sufficiently_decreased(phi_i, alpha_i)
             or ((phi_i >= phi_prev) and i > 0)):
            return results_for(
                    zoom(alpha_prev, alpha_i, phi_0, dphi_0, phi_prev,
                         phi_i, alpha_i
                    )
            )

        dphi_i = dphi_alpha(alpha_i)
        if curv_cond(dphi_i):
            return results_for(alpha_i)

        if dphi_i >= 0:
            return results_for(
                    zoom(alpha_i, alpha_prev, phi_0, dphi_0, phi_i,
                         phi_alpha_=phi_i, alpha_0_=alpha_i
                    )
            )
        prev_alpha = alpha
        alpha_i = min(fac * alpha_i, alpha_max)
    raise Exception("line_search() didn't converge in {i+1} cycles!")


def run_line_search(f, df, get_p, x0, alpha_0=1, alpha_max=1):
    x = x0
    alpha_prev = None
    gradients = list()
    step_dirs = list()

    for i in range(50):
        f_0 = f(x)
        grad = np.array(df(x))
        norm = np.linalg.norm(grad)
        print(f"{i:02d} norm(grad)={norm:.6f}")
        if norm < 1e-4:
            print("Converged!")
            break
        p = get_p(x)
        # p = -grad / np.linalg.norm(grad)

        gradients.append(grad)
        step_dirs.append(p)

        alpha_guess = 1
        if alpha_prev:
            # Try an improved guess for alpha
            # See [1], between Eq. (3.59) and (3.60), "Initial step length"
            numer = gradients[-2].dot(step_dirs[-2])
            denom = gradients[-1].dot(step_dirs[-1])
            alpha_guess = alpha_prev * numer/denom
            # Restrict the maximum value of the guessed alpha. Otherwise
            # huge alphas may result when the gradient reduction was big in
            # the previous cycle.
            alpha_guess = min(alpha_max, alpha_guess)
            print(f"\tusing alpha_guess={alpha_guess:.6f}")
            assert alpha_guess > 0

        # f_evals_prev = f_evals
        # df_evals_prev = df_evals
        alpha = line_search(f, df, x, p, f_0, df_0=grad, alpha=alpha_guess)
        # fc_ = f_evals - f_evals_prev
        # gc_ = df_evals - df_evals_prev
        # print(f"\talpha={alpha:.6f}, fc={fc_}, gc={gc_}")
        print(f"\talpha={alpha:.6f}")
        # sp_alpha, fc, gc, *_ = sp_line_search(f_, f_grad_, x, p)
        # if sp_alpha:
            # print(f"\tsp_alpha={sp_alpha:.6f}, fc={fc}, gc={gc}")
            # x = x + sp_alpha*p
        # else:
            # x = x + alpha*p
        x = x + alpha*p
        alpha_prev = alpha


def shift_hessian(hessian, beta=1e-3, max_cycles=10):
    w, v = np.linalg.eigh(hessian)
    print(f"\tw.min()={w.min():.4e}")
    if w.min() > 0:
        tau = 0
    else:
        tau = -(w.min()) + beta
    eye = np.eye(hessian.shape[0])
    for k in range(10):
        print(f"\t{k:02d}: tau={tau:.6e}")
        try:
            L = np.linalg.cholesky(hessian + tau*eye)
            return L @ L.T
        except np.linalg.LinAlgError:
            tau = max(2*tau, beta)
    raise Exception( "Shifting the hessian didn't succeed after "
                    f"{k+1} cycles.")


def rosenbrock():
    f_evals = 0
    df_evals = 0

    def f(x0, x1):
        nonlocal f_evals
        f_evals += 1
        return 100*(-x0**2 + x1)**2 + (1 - x0)**2

    def f_grad(x0, x1):
        nonlocal df_evals
        df_evals += 1
        return (-400*x0*(x1-x0**2) - 2 + 2*x0, 200*(x1-x0**2))

    def f_hess(x0, x1):
        return (
            (1200*x0**2 - 400*x1 + 2, -400*x0),
            (-400*x0, 200)
        )

    def f_(x):
        return f(*x)

    def f_grad_(x):
        return f_grad(*x)

    def get_p(x):
        grad = f_grad_(x)
        hess = np.array(f_hess(*x))
        newton_step = -np.linalg.inv(hess) @ grad
        return newton_step

    x0 = (1.2, 1.2)
    run_line_search(f_, f_grad_, get_p, x0)


def run():
    f_evals = 0
    df_evals = 0
    def f(x):
        nonlocal f_evals
        f_evals += 1
        return 2*x**4 + 5*x**3 - 2*x**2 + 10*x

    def f_grad(x):
        nonlocal df_evals
        df_evals += 1
        return (8*x**3 + 15*x**2 - 4*x + 10, )

    def f_(x):
        return f(*x)

    def f_grad_(x):
        return np.array(f_grad(*x))

    def get_p(x):
        grad = f_grad_(x)
        return -grad / np.linalg.norm(grad)

    x0 = (-3, )
    run_line_search(f_, f_grad_, get_p, x0)


def run_geom_opt():
    from pysisyphus.calculators.XTB import XTB
    from pysisyphus.helpers import geom_from_library

    geom = geom_from_library("birkholz/codeine.xyz")
    calc = XTB(pal=4)
    geom.set_calculator(calc)

    geom_opt = calc.run_opt(geom.atoms, geom.coords)
    geom_opt.set_calculator(XTB(pal=4))
    opt_hess = geom_opt.hessian
    wo, vo = np.linalg.eigh(opt_hess)
    # import pdb; pdb.set_trace()

    def f_(x):
        res = geom.get_energy_and_forces_at(x)
        return res["energy"]

    def f_grad_(x):
        res = geom.get_energy_and_forces_at(x)
        return -res["forces"]

    def get_p(x):
        grad = f_grad_(x)
        hess = geom.hessian
        hess_shifted = shift_hessian(hess)
        newton_step = -np.linalg.inv(hess) @ grad
        nn = np.linalg.norm(newton_step)
        newton_step_ = -np.linalg.inv(hess_shifted) @ grad
        nn_ = np.linalg.norm(newton_step_)
        # import pdb; pdb.set_trace()
        # return newton_step
        return newton_step_

    # run_line_search(f_, f_grad_, get_p, x0=geom.coords)
    for i in range(25):
        grad = geom.gradient
        norm = np.linalg.norm(grad)
        max_ = np.abs(grad).max()
        print(f"{i:02d} {norm:.4e} {max_:.4e}")
        hess = geom.hessian
        # for i in range(5):
            # hess_ = shift_hessian(hess, beta=10**(-i))
            # step = -np.linalg.pinv(hess_) @ grad
            # step_max = np.abs(step).max()
            # print(f"i={i}, step_max={step_max}")
        # import pdb; pdb.set_trace()
        # return
        hess_ = shift_hessian(hess, beta=1e-1)
        step = -np.linalg.pinv(hess_) @ grad
        new_coords = geom.coords + step
        geom.coords = new_coords
