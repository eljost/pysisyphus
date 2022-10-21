#!/usr/bin/env python

# [1] https://pubs.acs.org/doi/full/10.1021/acs.jctc.7b00788
#     libreta: Computerized Optimization and Code Synthesis for Electron
#     Repulsion Integral Evaluation
#     Jun Zhang, 2018
# [2] https://doi.org/10.1002/qua.560400604
#     Two-electron repulsion integrals over Gaussian s functions
#     Gill, Johnson, Pople, 1991

import argparse
import sys

import numpy as np
from scipy.special import factorial2

from pysisyphus.config import CONFIG_DIR, L_MAX, L_AUX_MAX


_BOYS_TABLE = np.load(
    CONFIG_DIR / "wavefunction/ints/boys_table_N_64_xasym_27.1_step_0.01.npy"
)
_BOYS_N_MAX = (
    max(L_MAX, L_AUX_MAX) * 3
)  # Currently, there are at most 3center integrals
_BOYS_X_ASYM = (
    27.0  # Beyond this x the Boys function can be calculated via an analytical formula
)
_BOYS_ADD_POINTS = 10  # Go a bit beyond X_ASYM
_BOYS_STEP = 0.01
_BOYS_FACTOR = 100  # int(1 / _BOYS_FACTOR)
_BOYS_XS = np.arange(0, _BOYS_X_ASYM + (_BOYS_ADD_POINTS + 1) * _BOYS_STEP, _BOYS_STEP)

PI = np.pi
FACT2 = np.array(
    [factorial2(N, exact=True) for N in range(2 * _BOYS_N_MAX)], dtype=np.int64
)


def make_boys_table():
    from scipy.integrate import quad

    def boys(t, x, N):
        return t ** (2 * N) * np.exp(-x * t**2)

    xs = _BOYS_XS
    boys_table = list()
    for N in range(_BOYS_N_MAX + 1):
        print(f"{N=}, calculating at {xs.size} points from {xs[0]} to {xs[-1]}")
        vals = list()
        max_abserr = 0.0
        for i, x in enumerate(xs):
            if i % 1000 == 0:
                print(f"\t{i} values")
            y, abserr = quad(boys, 0, 1, args=(x, N), epsabs=5e-13, epsrel=5e-13)
            max_abserr = max(abserr, max_abserr)
            vals.append(y)
        print(f"\tN={N}, max(abserr)={max_abserr:.8e}")
        boys_table.append(vals)
    boys_table = np.array(boys_table)
    fn = f"boys_table_N_{_BOYS_N_MAX}_xasym_{xs[-1]:.1f}_step_{_BOYS_STEP:.2f}.npy"
    np.save(fn, boys_table)
    return fn


def neville(x, xs_table, table, step, points=5, factor=None):
    """Slow, recursive implementation of Neville interpolation.

    We multiply 'x' by int(1/step), so we only have to deal with integer arithmetic, to
    determine the first entry from the table. 'x_closest' is always chosen in a way, that
    'x' is contained in the 'xs' interval. Just to be sure, one can also suppy the 'factor'
    to this method."""

    if factor is None:
        factor = int(1 / step)

    # Determine closest x
    x_closest = int(x * factor)
    indices = np.arange(points) + x_closest
    xs = xs_table[indices]
    ys = table[indices]
    p_ = {(i, i): ys[i] for i in range(points)}

    def p(i, j):
        key = (i, j)
        if key not in p_:
            xi, xj = xs[[i, j]]
            p_[key] = ((x - xj) * p(i, j - 1) - (x - xi) * p(i + 1, j)) / (xi - xj)
        return p_[key]

    return p(0, points - 1)


def neville_gen5(x, xs_table, table, step, factor=None):
    """Code-generated Neville interpolation using 5 points."""

    points = 5
    if factor is None:
        factor = int(1 / step)

    # Determine closest x
    x_closest = int(x * factor)
    indices = np.arange(points) + x_closest
    xs = xs_table[indices]
    ys = table[indices]

    x0 = -xs[4]
    x1 = x - xs[0]
    x2 = -xs[1]
    x3 = x + x2
    x4 = -xs[2]
    x5 = x + x4
    x6 = -xs[3]
    x7 = x + x6
    x8 = x + x0
    x9 = (x5 * ys[3] - x7 * ys[2]) / (-x6 - xs[2])
    x10 = (x3 * ys[2] - x5 * ys[1]) / (-x4 - xs[1])
    x11 = (-x10 * x7 + x3 * x9) / (-x6 - xs[1])

    return (
        x1
        * (
            -x11 * x8
            + x3
            * (x5 * (x7 * ys[4] - x8 * ys[3]) / (-x0 - xs[3]) - x8 * x9)
            / (-x0 - xs[2])
        )
        / (-x0 - xs[1])
        - x8
        * (
            x1 * x11
            - x7
            * (x1 * x10 - x5 * (x1 * ys[1] - x3 * ys[0]) / (-x2 - xs[0]))
            / (-x4 - xs[0])
        )
        / (-x6 - xs[0])
    ) / (-x0 - xs[0])


def neville_boys(N, x):
    """Wrapper for Boys function from Neville-interpolation."""
    return neville_gen5(
        x, _BOYS_XS, _BOYS_TABLE[N], step=_BOYS_STEP, factor=_BOYS_FACTOR
    )


def factorial_boys(N, x):
    """Boys function for (big) x > 27.0 as described in the SI of [1]. See also [2].
                  ___________
     -N - 1      ╱  -2⋅N - 1
    2      ⋅√π⋅╲╱  x         ⋅(2⋅N - 1)!

    """
    return FACT2[max(0, 2 * N - 1)] / 2 ** (N + 1) * np.sqrt(PI / x ** (2 * N + 1))


def boys(N, xs):
    """Wrapper for Boys function calculation.

    Supports scalar values and np.ndarray for 'xs'.
    """

    def func(N, x):
        return neville_boys(N, x) if (x <= _BOYS_X_ASYM) else factorial_boys(N, x)

    if isinstance(xs, np.ndarray):
        boys_list = list()

        with np.nditer(xs) as it:
            for x in it:
                boys_list.append(func(N, x))
        boys_table = np.reshape(boys_list, xs.shape)
        return boys_table
    else:
        return func(N, xs)


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--table", action="store_true")
    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    if args.table:
        fn = make_boys_table()
        print(f"Wrote Boys table to '{fn}'.")


if __name__ == "__main__":
    run()
