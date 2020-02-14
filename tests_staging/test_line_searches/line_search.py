#!/usr/bin/env python3

# Based on
# [1] Nocedal - Numerical Optimization, chapter 3. Line Search Methods

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import line_search as sp_line_search

from pysisyphus.optimizers.line_searches import wolfe


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
        # import pdb; pdb.set_trace()

        kwargs = {
            "f": f,
            "df": df,
            "x0": x,
            "p": p,
            "f0": f_0,
            "g0": grad,
            "alpha_init": alpha_guess,
        }
        alpha, f_, df_ = wolfe(**kwargs)

        # alpha, f_, df_ = wolfe(f, df, x, p, f_0, df_0=grad, alpha=alpha_guess)
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


if __name__ == "__main__":
    run()
