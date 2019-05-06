#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root_scalar

from pysisyphus.helpers import geom_from_library, shake_coords
from pysisyphus.calculators.XTB import XTB


np.set_printoptions(suppress=True, precision=4, linewidth=120)


def run():
    # geom = geom_from_library("benzene.xyz")
    geom = geom_from_library("codein.xyz")
    geom.coords = shake_coords(geom.coords, seed=25032019)
    calc = XTB(charge=0, mult=1, pal=2)
    geom.set_calculator(calc)

    # g = geom.gradient
    # np.savetxt("gradient", g)
    # H = geom.hessian
    # np.savetxt("hessian", H) 
    H = np.loadtxt("hessian")
    g = np.loadtxt("gradient")
    trust = 0.3


    newton_step = -np.linalg.pinv(H).dot(g)
    newton_norm = np.linalg.norm(newton_step)
    print("newton step norm", newton_norm)
    vals, vecsT = np.linalg.eigh(H)
    neg_inds = vals < 0
    H_ = np.diag(vals)
    g_ = vecsT.T.dot(g)
    f1 = g_[0]
    print(f"f1 is {f1}")

    def step_(lambda_):
        _ = g_/(vals+lambda_)
        return _ * vecsT
    step_(0.5)

    def step_norm(lambda_):
        return np.linalg.norm(step_(lambda_))

    def on_sphere(lambda_):
        return step_norm(lambda_) - trust

    def on_sphere_linear(lambda_):
        return 1/trust - 1/step_norm(lambda_)

    _b1 = -vals[0]
    x0 = _b1 + 1e-3
    x1 = x0 + 1e-3

    # bracket = [_b1, float("inf")] doesn't work, use big number instead
    upper_bracket = 1e10
    # Hessian  is positive definite
    if neg_inds.sum() == 0:
        bracket = [0, upper_bracket]
    # Hessian has negative eigenvalues, is not positive definitie
    else:
        bracket = [_b1+1e-6, upper_bracket]
    # sol = root_scalar(on_sphere, x0=x0, x1=x1, bracket=bracket)
    sol_lin = root_scalar(on_sphere_linear, x0=x0, x1=x1, bracket=bracket)
    print(sol_lin)

    import pdb; pdb.set_trace()
    # g_[0] = 0.0001
    # g_[0] = 0.0005
    # Case 3 if g_[0] <= 1e-4
    lambdas = np.linspace(0, 1, 500)
    step_norms = list()
    eye = np.eye(H_.shape[0])
    for l in lambdas:
        step = -np.linalg.pinv(H_ + l*eye).dot(g_)
        norm = np.linalg.norm(step)
        step_norms.append(norm)
        # print(f"l={l:.4e}, norm={norm:.6f}")
    fig, ax = plt.subplots()
    ax.plot(lambdas, step_norms)
    for neg_eigval in vals[neg_inds]:
        ax.axvline(-neg_eigval, linestyle="--", color="red")
    ax.axhline(trust, linestyle=":", color="k")
    ax.set_ylim(0, 5)

    plt.show()


if __name__ == "__main__":
    run()
