import functools
from math import prod

import numpy as np
from numpy.typing import NDArray
from scipy.special import factorial2

from pysisyphus.wavefunction.helpers import canonical_order


# @functools.cache
def get_lmn_factors(L: int):
    lmns = canonical_order(L)
    lmn_factors = np.zeros(len(lmns))
    for i, lmn in enumerate(lmns):
        lmn_factors[i] = prod([factorial2(2 * am - 1) for am in lmn])
    lmn_factors = 1 / np.sqrt(lmn_factors)
    return lmn_factors


def norm_cgto_lmn(coeffs: NDArray[float], exps: NDArray[float], L: int):
    N = 0.0
    for i, expi in enumerate(exps):
        for j, expj in enumerate(exps):
            tmp = coeffs[i] * coeffs[j] / (expi + expj) ** (L + 1.5)
            tmp *= np.sqrt(expi * expj) ** (L + 1.5)
            N += tmp
    N = np.sqrt(exps ** (L + 1.5) / (np.pi**1.5 / 2**L * N))
    # Or w/ array broadccasting
    # N = coeffs[None, :] * coeffs[:, None] / (exps[None, :] + exps[:, None]) ** (L + 1.5)
    # N = (N * np.sqrt((exps[None, :] * exps[:, None]) ** (L + 1.5))).sum()
    # N = np.sqrt(exps ** (L + 1.5) / (np.pi**1.5 / 2**L * N))

    mod_coeffs = N * coeffs
    lmn_factors = get_lmn_factors(L)

    return mod_coeffs, lmn_factors
