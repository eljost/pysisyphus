#!/usr/bin/env python3
# Johannes Steinmetzer, April 2019

# See [1] https://pubs.acs.org/doi/pdf/10.1021/j100247a015
#         Banerjee, 1985
#     [2]
#

import matplotlib.pyplot as plt
import numpy as np
import sympy as sym


def make_funcs():
    x, y, = sym.symbols("x y")
    f_ = (1 - y)*x**2*sym.exp(-x**2) + (1/2)*y**2
    f = sym.lambdify((x, y), f_)
    g_ = sym.derive_by_array(f_, (x, y))
    g = sym.lambdify((x, y), g_)
    H_ = sym.derive_by_array(g_, (x, y))
    H = sym.lambdify((x, y), H_)
    return f, g, H


def plot(f, g, H, xs, ys):
    X, Y = np.meshgrid(xs, ys)
    Z = f(X, Y)
    levels = np.linspace(0, 2, 75)
    # fig, ax = plt.subplots(figsize=(12, 8))
    # cf = ax.contour(X, Y, Z, levels=levels)
    # fig.colorbar(cf)
    # plt.show()
    
    neg_eigvals = list()
    grads = list()
    for x_ in xs:
        for y_ in ys:
            hess = H(x_, y_)
            eigvals = np.linalg.eigvals(hess)
            if eigvals.min() < 0:
                neg_eigvals.append((x_, y_))
            grad = np.linalg.norm(g(x_, y_))
            grads.append(grad)
    neg_eigvals = np.array(neg_eigvals)
    grads = np.array(grads)
    
    fig, ax = plt.subplots(figsize=(12, 8))
    cf = ax.contour(X, Y, Z, levels=levels)
    ax.scatter(*neg_eigvals.T, c="r", s=15, label="neg. eigval.")
    ax.scatter(X.T, Y.T, c="b", s=5*grads, label="norm(grad)")
    ax.legend()
    fig.colorbar(cf)
    plt.show()


def prfo(x, H_getter, grad_getter):

    fg = lambda x: -np.array(grad_getter(*x))
    Hg = lambda x:  np.array(H_getter(*x))
    f = fg(x)
    H = Hg(x)
    eigvals, eigvecs = np.linalg.eigh(H)
    neg_eigvals = eigvals < 0
    assert neg_eigvals.sum() >= 1
    print(f"found {neg_eigvals.sum()} negative eigenvalues")

    # Transform to eigensystem of hessian
    f_trans = eigvecs.T.dot(f)
    
    mu = 0
    max_mat = np.array(((eigvals[mu], -f_trans[mu]),
                       (-f_trans[mu], 0)))
    min_mat = np.bmat((
        (np.diag(eigvals[1:]), -f_trans[1:,None]),
        (-f_trans[None,1:], [[0]])
    ))
    
    # Scale eigenvectors of the largest (smallest) eigenvector 
    # of max_mat (min_mat) so the last item is 1.
    max_evals, max_evecs = np.linalg.eigh(max_mat)
    # Eigenvalues and -values are sorted, so we just use the last
    # eigenvector corresponding to the biggest eigenvalue.
    max_step = max_evecs.T[-1]
    lambda_max = max_step[-1]
    max_step = max_step[:-1] / lambda_max

    min_evals, min_evecs = np.linalg.eigh(min_mat)
    # Again, as everything is sorted we use the (smalelst) first eigenvalue.
    min_step = np.asarray(min_evecs.T[0]).flatten()
    lambda_min = min_step[-1]
    min_step = min_step[:-1] / lambda_min

    # Still in the hessian eigensystem
    prfo_step = np.zeros_like(f)
    prfo_step[0] = max_step[0]
    prfo_step[1:] = min_step
    # Transform back
    step = eigvecs.dot(prfo_step)
    norm = np.linalg.norm(step)
    if norm > 0.1:
        step = 0.1 * step / norm
    return step


def run():
    x0 = (0.5, 0.2)
    f, g, H = make_funcs()
    xs = np.linspace(-1.3, 1.3, 50)
    ys = np.linspace(-0.7, 1.9, 50)
    plot(f, g, H, xs, ys)
    prfo(x0, H, g)
    fq = -np.array(g(*x0))
    hess = np.array(H(*x0))

    x = x0
    opt_xs = [x, ]

    for i in range(15):
        step = prfo(x, H, g)
        print("norm(step)", np.linalg.norm(step))
        grad = g(*x)
        gn = np.linalg.norm(grad)
        if gn < 1e-5:
            print("Converged")
            break
        x_new = x + step
        opt_xs.append(x_new)
        x = x_new
    opt_xs = np.array(opt_xs)

    X, Y = np.meshgrid(xs, ys)
    Z = f(X, Y)
    levels = np.linspace(0, 2, 75)
    fig, ax = plt.subplots()
    cf = ax.contour(X, Y, Z, levels=levels)
    ax.plot(*opt_xs.T, "ro-", label="TSopt")
    fig.colorbar(cf)
    ax.set_xlim(xs.min(), xs.max())
    ax.set_ylim(ys.min(), ys.max())
    plt.show()

    
if __name__ == "__main__":
    run()