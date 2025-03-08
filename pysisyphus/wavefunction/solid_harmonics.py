from math import floor

import numpy as np
import scipy as sp


"""
Regular solid harmonics, based on formulas from wikipedia
    https://en.wikipedia.org/wiki/Solid_harmonics#Real_form
"""

# Imaginary unit
Im = 1.0j
# Factorial function
fact = sp.special.factorial
# Binomical coefficient
binom = sp.special.binom
HALF = 0.5

SQRT3 = np.sqrt(3.0)
HALF_SQRT3 = SQRT3 / 2.0


def Am(m, x, y):
    return HALF * ((x + Im * y) ** m + (x - Im * y) ** m)


def Bm(m, x, y):
    return 1 / (2 * Im) * ((x + Im * y) ** m - (x - Im * y) ** m)


def gamma(l, m, k):
    return (
        (-1) ** k
        * 2 ** (-l)
        * binom(l, k)
        * binom(2 * l - 2 * k, l)
        * fact(l - 2 * k)
        / fact(l - 2 * k - m)
    )


def zpart(l, m, x, y, z):
    limit = floor((l - m) / 2)
    sum_ = 0.0
    r = np.sqrt(x**2 + y**2 + z**2)
    for k in range(limit + 1):
        sum_ += gamma(l, m, k) * r ** (2 * k) * z ** (l - 2 * k - m)
    return sum_


def C(l, m, x, y, z):
    krond = int(m == 0)
    prefact = ((2 - krond) * fact(l - m) / fact(l + m)) ** HALF
    result = prefact * zpart(l, m, x, y, z) * Am(m, x, y)
    return result


def S(l, m, x, y, z):
    prefact = (2 * fact(l - m) / fact(l + m)) ** HALF
    result = prefact * zpart(l, m, x, y, z) * Bm(m, x, y)
    return result


def Rlm(l, m, x, y, z):
    """Real regular solid harmonics w/ Cartesian arguments."""
    assert abs(m) <= l, f"{m=} must be in the range {-l=} <= m <= {l=}!"
    # Some hardcoded shortcuts below
    if l == 0:
        return 1.0
    elif l == 1:
        # 0, 1, -1
        return (z, x, y)[m]
    elif l == 2:
        # 0, 1, 2, -2, -1
        return (
            lambda x, y, z: 0.5 * (3 * z**2 - x**2 - y**2 - z**2),  # 0
            lambda x, y, z: SQRT3 * x * z,  # 1
            lambda x, y, z: HALF_SQRT3 * (x**2 - y**2),  # 2
            lambda x, y, z: SQRT3 * x * y,  # -2
            lambda x, y, z: SQRT3 * y * z,  # -1
        )[m](x, y, z)

    if m < 0:
        func = S
        m = abs(m)
    else:
        func = C
    return np.real(func(l, m, x, y, z))
