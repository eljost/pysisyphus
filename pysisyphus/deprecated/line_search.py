#!/usr/bin/env python3

import numpy as np

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


if __name__ == "__main__":
    run()
