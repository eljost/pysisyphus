#!/usr/bin/env python3

from math import sqrt

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import root_scalar

from pysisyphus.optimizers.hessian_updates import bfgs_update
from pysisyphus.helpers import geom_from_library, shake_coords, geom_from_xyz_file
from pysisyphus.calculators.XTB import XTB


np.set_printoptions(suppress=True, precision=4, linewidth=120)


# TODO: determine lambda for the 3 cases and calculate the step
# at the end, not in three places.
def rsa(H, g, trust=0.3):
    vals, vecsT = np.linalg.eigh(H)
    neg_inds = vals < 0
    neg = neg_inds.sum()

    H_ = np.diag(vals)
    g_ = vecsT.T.dot(g)

    def get_step(lambda_):
        """Returns step for a given lambda"""
        # _ = -g_/(vals+lambda_)
        _ = -g_/(vals+1e-12+lambda_)
        return vecsT.dot(_)

    if neg == 0:
        # print("Hessian is positive definite.")
        # print("Checking pure QN-step ... ", end="")
        step = get_step(0)
        step_norm = np.linalg.norm(step)
        if step_norm <= trust:
            # print("is fine")
            return step
        # print(f"norm(QN-step)={step_norm:.6f} is too big.")
    else:
        # Hessian is not positive definite
        # print("Hessian is not positive definite.")
        smallest_eigval = vals[0]
        print(f"Smallest eigenvalue is {smallest_eigval:.6f}")

    def on_sphere_linear(lambda_):
        return 1/trust - 1/np.linalg.norm(get_step(lambda_))

    _b1 = -vals[0]
    x0 = _b1 + 1e-3
    x1 = x0 + 1e-3

    # Defining a bracket using infinity (float("inf")) doesn't seem to work.
    # Instead we use a big number.
    upper_bracket = 1e10
    # Hessian  is positive definite
    if neg == 0:
        bracket = [0, upper_bracket]
    # Hessian has negative eigenvalues, is not positive definitie
    else:
        bracket = [_b1+1e-6, upper_bracket]
    sol = root_scalar(on_sphere_linear, x0=x0, x1=x1, bracket=bracket)
    if not sol.converged:
        raise Exception("Root search did not converge!")
    lambda_ = sol.root
    # Check shifted hessian for positive definiteness by adding the shift
    # to the smallest eigenvalue and see if this is > 0. If so we can use this
    # step (Case II).
    if vals[0] + lambda_ > 0:
        step = get_step(lambda_)
        step_norm = np.linalg.norm(step)
        print(f"Found valid step with λ={lambda_:.6f} and norm={step_norm:.6f}")
        return step

    import pdb; pdb.set_trace()
    print(f"Shifted hessian (λ={lambda_:.6f} is not positive definite!")
    print("Determining new step using second parameter τ.")
    # Shifted hessian is not positive definite (Case III).
    lambda_ = vals[0]
    # frac_sum = np.sum(g_[1:] / (vals[1:] - vals[0]))
    frac_sum = np.sum(g_[1:] / (vals[1:] - lambda_))
    tau = sqrt(trust**2 - frac_sum**2)
    print(f"τ={tau:.6f}")

    step = get_step(lambda_) + tau*g_[0]
    return step


def run():
    # pysis_in = np.loadtxt("ref_rsa/in_hess")
    # this_in = np.loadtxt("test_in_hess")
    # geom = geom_from_library("codein.xyz")
    # geom = geom_from_xyz_file("ref_rsa/codeine.xyz")
    geom = geom_from_library("birkholz/vitamin_c.xyz")
    # geom.coords = shake_coords(geom.coords, seed=25032019)
    calc = XTB(charge=0, mult=1, pal=4)
    geom.set_calculator(calc)
    from pysisyphus.optimizers.RSAlgorithm import RSAlgorithm
    rsa_kwargs = {
        "hessian_recalc": 5,
    }
    opt = RSAlgorithm(geom, **rsa_kwargs)
    opt.run()
    return
    # g = geom.gradient
    # np.savetxt("gradient", g)
    # H = geom.hessian
    H = np.eye(geom.coords.size)
    np.savetxt("test_in_hess", H) 
    # H = np.loadtxt("hessian")
    # g = np.loadtxt("gradient")
    # H = geom.hessian
    # from pysisyphus.calculators.AnaPot import AnaPot
    # pot = AnaPot()
    # geom = AnaPot.get_geom((0, 3, 0))
    # H = geom.hessian
    # H = np.eye(geom.coords.shape[0]) / 2
    trust = 0.3
    max_cycles = 50

    coords = list()
    steps = list()
    grads = list()
    # import pdb; pdb.set_trace()
    for i in range(max_cycles):
        coords.append(geom.coords.copy())
        g = geom.gradient
        grads.append(g)
        grad_norm = np.linalg.norm(g)
        if i > 0:
            dg = grads[-1] - grads[-2]
            dx = steps[-1]
            dH, _ = bfgs_update(H, dx, dg)
            H = H + dH
            # H = geom.hessian
        rms_g = np.sqrt(np.mean(g**2))
        print(f"Cycle {i:02d}: norm(g)={grad_norm:.4e}, rms(g)={rms_g:.6f}")
        if grad_norm < 1e-3:
            print("Converged")
            break
        # import pdb; pdb.set_trace()
        step = rsa(H, g, trust)
        new_coords = geom.coords + step
        geom.coords = new_coords
        steps.append(step)
    with open("opt.xyz", "w") as handle:
        handle.write(geom.as_xyz())
    coords = np.array(coords)
    # pot.plot()
    # ax = pot.ax
    # ax.plot(*coords.T[:2], "bo-")
    plt.show()


if __name__ == "__main__":
    run()
