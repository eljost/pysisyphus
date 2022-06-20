# [1] https://pubs.acs.org/doi/full/10.1021/acs.jctc.7b00788
#     libreta: Computerized Optimization and Code Synthesis for Electron
#     Repulsion Integral Evaluation
#     Jun Zhang, 2018

from math import ceil, log
import numpy as np
from scipy.special import factorial2


from pysisyphus.config import CONFIG_DIR


_BOYS_TABLE = np.load(CONFIG_DIR / "wavefunction/boys_table.npy")
_BOYS_STEP = 0.01
_BOYS_NDIGITS = 2
# _BOYS_N_MAX = _BOYS_TABLE.shape[0] - 1
PI = np.pi


def neville(x, table, step, points=5, ndigits=None):
    # Determine significant places
    if ndigits is None:
        ndigits = ceil(abs(log(step, 10)))

    # Determine closest x
    res = abs(x) % step
    sign = 1 if (res / step < 0.5) else -1
    x_closest = round(x + sign * res, ndigits=ndigits)

    xs = x_closest + step * np.arange(points)
    indices = (np.round(xs, ndigits) / step).astype(int)
    ys = table[indices]
    p_ = {(i, i): ys[i] for i in range(points)}

    def p(i, j):
        key = (i, j)
        if key not in p_:
            xi, xj = xs[[i, j]]
            p_[key] = ((x - xj) * p(i, j - 1) - (x - xi) * p(i + 1, j)) / (xi - xj)
        return p_[key]

    return p(0, points - 1)


def neville_boys(N, x):
    return neville(x, _BOYS_TABLE[N], step=_BOYS_STEP, ndigits=_BOYS_NDIGITS)


def factorial_boys(N, x):
    return factorial2(2 * N - 1) / 2 ** (N + 1) * np.sqrt(PI / x ** (2 * N + 1))


def boys(N, x):
    result = neville_boys(N, x) if x <= 27.0 else factorial_boys(N, x)
    return result
