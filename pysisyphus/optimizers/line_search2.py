#!/usr/bin/env python3

from collections import namedtuple
import logging
from math import sqrt
from pprint import pprint

import numpy as np
import sympy as sym


logger = logging.getLogger("optimizer")


def log(msg):
    logger.debug(msg)


def gen_solutions():
    """
    Given two energies (e0, e1) and corresponding gradients (g0, g1) we
    can (try to) fit a quartic polynomial
        f(x) = a0 + a1*x + a2*x**2 + a3*x**3 + a4*x**4
    s.t. the constraint f''(x) >= 0, with the equality being fullfilled
    at only one point.
    There are five unknowns (a0 - a4) to be determined. Four equations can
    be derived from f(x) and its first derivative
        f'(x) = a1 + 2*a2*x + 3*a3*x**2 + 4*a4*x**3 .
    With (e0, g0) being given at x=0 and (e1, g1) being given at x=1 we can
    setup the following equations:
        f (0) = a0                          (1)
        f'(0) = a1                          (2)
    using e0 and g0 at x=0, and
        f (1) = a0 + a1 + a2 + a3 + a4      (3)
        f'(1) = a1 + 2*a2 + 3*a3 + 4*a4 .   (4)
    The missing last equation can be derived from the constraint. The second
    derivative of f(x) is
        f''(x) = 2*a2 + 6*a3*x + 12*a4*x**2
    and shall be positive except at one point where it is allowed to be 0, that
    its two roots (f''(x) = 0) must be degenerate. This is fullfilled when the
    discriminant D of the quadratic polynomial a*x**2 + b*x + c is zero.
        D = b**2 – 4*a*c = 0
    With
        a = 12*a4
        b =  6*a3
        c =  2*a2
    we get
        0 = (6*a3)**2 - 4*12*a4*2*a2
        0 = 36*a3**2 - 96*a4*a2
        0 = 3*a3**2 - 8*a4*a2               (5)
        or
        a4 = 3/8 * a3**2 / a2
    Using (1) - (5) we can solve the set of equations for a0 - a4.
    """

    e0, e1, g0, g1, a0, a1, a2, a3 = sym.symbols("e0 e1 g0 g1 a:4")

    a4 = sym.Rational(3, 8) * a3**2 / a2
    s0, s1 = sym.solve((e0-a0,
                        g0-a1,
                        e1-a0-a1-a2-a3-a4,
                        g1-a1-2*a2-3*a3-4*a4,
                        3*a3**2 - 8*a2*a4),
                        (a0, a1, a2, a3)
    )
    print("Solution 0")
    print("\t", s0)
    print()
    print("Solution 1")
    print("\t", s1)
    print()
    # There will be two solutions (s0, s1), both containing a big sqrt(...) term
    # that can be computed once and reused.
    s0_cse = sym.cse(s0)
    s1_cse = sym.cse(s1)
    print("Solution 0 after CSE")
    pprint(s0_cse)
    print("Solution 1 after CSE")
    pprint(s1_cse)
    print()
    # The terms in the sqrt-term correspond to binomial expansions and can be further
    # simplified.
    ref_term = -12*e0**2 + 24*e0*e1 - 12*e0*g0 - 12*e0*g1 - 12*e1**2 + \
                12*e1*g0 + 12*e1*g1 - 2*g0**2 - 8*g0*g1 - 2*g1**2
    sqrt_term = -2*(6*(e0-e1)**2 + 6*(e0-e1)*(g0+g1) + (g0+g1)**2 + 2*g0*g1)
    assert sym.simplify(sym.expand(sqrt_term) - ref_term) == 0


def get_minimum(poly):
    roots = np.roots(np.polyder(poly))
    real_roots = np.real(roots[np.isreal(roots)])
    vals = poly(real_roots)
    min_ind = vals.argmin()
    min_root = real_roots[vals.argmin()]
    min_val = vals[min_ind]
    return min_root, min_val


FitResult = namedtuple("FitResult", "x y polys")


def quintic_fit(e0, e1, g0, g1, H0, H1):
    a = -H0/2 + H1/2 - 6*e0 + 6*e1 - 3*g0 - 3*g1
    b = 3*H0/2 - H1 + 15*e0 - 15*e1 + 8*g0 + 7*g1
    c = -3*H0/2 + H1/2 - 10*e0 + 10*e1 - 6*g0 - 4*g1
    d = H0/2
    e = g0
    f = e0

    poly = np.poly1d((a, b, c, d, e, f))
    try:
        mr, mv = get_minimum(poly)
    except ValueError:
        return None

    fit_result = FitResult(mr, mv, (poly, ))
    return fit_result


def quartic_fit(e0, e1, g0, g1):
    """See gen_solutions() for derivation."""
    a0 = e0
    a1 = g0
    try:
        sqrt_term = sqrt(-2*(6*(e0-e1)**2 + 6*(e0-e1)*(g0+g1) + (g0+g1)**2 + 2*g0*g1))
    except ValueError:
        # In these cases there is no intermediate minimum between 0 and 1 and the term
        # under the square root becomes negative.
        return None

    a2_pre = -3*(e0 - e1) - 5*g0/2 - g1/2
    a3_pre = 2*e0 - 2*e1 + 2*g0

    def get_poly(a3, a2, a1, a0):
        a4 = 3/8 * a3**2 / a2
        return np.poly1d((a4, a3, a2, a1, a0))

    a2 = a2_pre - sqrt_term/2
    a3 = a3_pre + sqrt_term
    poly0 = get_poly(a3, a2, a1, a0)

    a2 = a2_pre + sqrt_term/2
    a3 = a3_pre - sqrt_term
    poly1 = get_poly(a3, a2, a1, a0)


    mr0, mv0 = get_minimum(poly0)
    mr1, mv1 = get_minimum(poly1)

    mr, mv = (mr0, mv0) if mv0 < mv1 else (mr1, mv1)

    # Shorter sympy implementation. Probably slower? But shouldn't matter...
    # a0, a1, a2, a3 = sym.symbols("a:4")
    # a4 = sym.Rational(3, 8) * a3**2 / a2
    # s0, s1 = sym.solve((e0-a0,
                        # g0-a1,
                        # e1-a0-a1-a2-a3-a4,
                        # g1-a1-2*a2-3*a3-4*a4,
                        # 3*a3**2 - 8*a2*a4),
                        # (a0, a1, a2, a3)
    # )
    # N = lambda exprs: [sym.N(expr) for expr in exprs]
    # sym_poly0 = get_poly(*N(s0[::-1]))
    # sym_poly1 = get_poly(*N(s1[::-1]))

    fit_result = FitResult(mr, mv, (poly0, poly1))
    return fit_result


def cubic_fit(e0, e1, g0, g1):
    # # Shorter sympy implementation. Probably slower? But shouldn't matter...
    # # Ok it is really slow ... and it's gone.
    # a0, a1, a2, a3 = sym.symbols("a:4")
    # s = sym.solve((e0-a0,
                   # g0-a1,
                   # e1-a0-a1-a2-a3,
                   # g1-a1-2*a2-3*a3),
                   # (a0, a1, a2, a3),
    # )
    # coeffs = [float(sym.N(expr)) for expr in (s[a3], s[a2], s[a1], s[a0])]
    d = e0
    c = g0
    b = -(g1 + 2*g0 + 3*e0 - 3*e1)
    a = 2*(e0 - e1) + g0 + g1
    # np.testing.assert_allclose([a, b, c, d], coeffs, atol=1e-10)
    poly = np.poly1d((a, b, c, d))
    try:
        mr, mv = get_minimum(poly)
    except ValueError:
        return None

    fit_result = FitResult(mr, mv, (poly, ))
    return fit_result


def poly_line_search(cur_energy, prev_energy, cur_grad, prev_grad, prev_step, coords):
    # TODO: always call line_search? Probably, because we can also extrapolate
    # in a linesearch.

    # energy_increased = (cur_energy - prev_energy) > 0.
    # if not energy_increased:
        # return cur_grad

    # Generate directional gradients by projecting them on the previous step.
    prev_grad_proj = prev_step @ prev_grad
    cur_grad_proj =  prev_step @ cur_grad
    cubic_result = cubic_fit(prev_energy, cur_energy, prev_grad_proj, cur_grad_proj)
    quartic_result = quartic_fit(prev_energy, cur_energy, prev_grad_proj, cur_grad_proj)
    # TODO: add quintic, but then we would have to save the hessians.

    accept = {
        # They way cubic is defined now it is never accepted and this is
        # probably better, because it doesn't seem to improve the optimization.
        "cubic": lambda x: (x > 2) and (x < 1),  # lgtm [py/redundant-comparison]
        "quartic": lambda x: (x > 0.) and (x <= 2),
    }
    fit_result = None
    if quartic_result and accept["quartic"](quartic_result.x):
        fit_result = quartic_result
        deg = "quartic"
    elif cubic_result and accept["cubic"](cubic_result.x):
        fit_result = cubic_result
        deg = "cubic"
    # else:
        # Midpoint fallback as described by gaussian?

    fit_step = None
    fit_grad = None
    fit_energy = None
    if fit_result and fit_result.y < prev_energy:
        x = fit_result.x
        fit_energy = fit_result.y
        log(f"Did {deg} interpolation with x={x:.6f}.")
        # Interpolate step and gradient
        fit_step = (1-x) * -prev_step
        fit_grad = (1-x)*prev_grad + x*cur_grad
    return fit_step, fit_grad, fit_energy
