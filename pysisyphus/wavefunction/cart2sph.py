#!/usr/bin/env python

# [1] https://doi.org/10.1002/qua.560540202
#     Transformation between Cartesian and pure spherical harmonic Gaussians
#     Schlegel, Frisch, 1995

import argparse
import functools
from cmath import sqrt
from math import floor
import sys
from typing import Dict, List, Tuple

import mpmath
from mpmath import factorial as mp_fact, sqrt as mp_sqrt, binomial as mp_binom
import numpy as np
from scipy.special import factorial as fact
from scipy.special import binom

from pysisyphus.wavefunction.helpers import canonical_order
from pysisyphus.wavefunction.normalization import get_lmn_factors


ZERO_THRESH = 1e-14


@functools.cache
def cart2sph_coeff(lx: int, ly: int, lz: int, m: int) -> complex:
    """Coefficient to convert Cartesian to spherical harmonics.

    Eq. (15) in [1].

    Paramters
    ---------
    lx
        Cartesian quantum number.
    ly
        Cartesian quantum number.
    lz
        Cartesian quantum number.
    m
        Magnetic quantum number.

    Returns
    -------
    coeff
        Contribution of the given Cartesian basis function to the spherical
        harmonic Y^(lx + ly + lz)_m.
    """

    l = lx + ly + lz
    assert -l <= m <= l

    j = (lx + ly - abs(m)) / 2
    if (j % 1) != 0:  # Return when j is not integer
        return 0
    j = int(j)

    lfact = fact(l)  # l!
    lmm = l - abs(m)  # (l - |m|)
    sign = 1 if m >= 0 else -1

    coeff = 0.0j
    for i in range(floor(lmm / 2) + 1):
        i_fact = (
            binom(l, i)
            * binom(i, j)
            * (-1) ** i
            * fact(2 * l - 2 * i)
            / fact(lmm - 2 * i)
        )
        k_sum = 0.0
        for k in range(j + 1):
            prefact = (-1) ** (sign * (abs(m) - lx + 2 * k) / 2)
            k_sum += prefact * binom(j, k) * binom(abs(m), lx - 2 * k)
        coeff += i_fact * k_sum
    coeff *= sqrt(
        fact(2 * lx)
        * fact(2 * ly)
        * fact(2 * lz)
        * lfact
        * fact(lmm)
        / (fact(2 * l) * fact(lx) * fact(ly) * fact(lz) * fact(l + abs(m)))
    )
    coeff *= 1 / (2**l * lfact)

    # print(f"{lx=}, {ly=}, {lz=}, {m=}, {coeff}")
    return coeff


@functools.cache
@mpmath.workdps(24)
def mp_cart2sph_coeff(lx: int, ly: int, lz: int, m: int) -> mpmath.mpc:
    """Higher precision Cartesian to Spherical transformation coefficients.

    Eq. (15) in [1].

    Paramters
    ---------
    lx
        Cartesian quantum number.
    ly
        Cartesian quantum number.
    lz
        Cartesian quantum number.
    m
        Magnetic quantum number.

    Returns
    -------
    coeff
        Contribution of the given Cartesian basis function to the spherical
        harmonic Y^(lx + ly + lz)_m.
    """

    l = lx + ly + lz
    assert -l <= m <= l

    j = (lx + ly - abs(m)) / 2
    if (j % 1) != 0:  # Return when j is not integer
        return mpmath.mp.mpc(0)
    j = int(j)

    lfact = mp_fact(l)  # l!
    lmm = l - abs(m)  # (l - |m|)
    sign = 1 if m >= 0 else -1

    coeff = mpmath.mp.mpc(0.0j)
    for i in range(floor(lmm / 2) + 1):
        i_fact = (
            mp_binom(l, i)
            * mp_binom(i, j)
            * (-1) ** i
            * mp_fact(2 * l - 2 * i)
            / mp_fact(lmm - 2 * i)
        )
        k_sum = 0.0
        for k in range(j + 1):
            prefact = (-1) ** (sign * (abs(m) - lx + 2 * k) / 2)
            k_sum += prefact * mp_binom(j, k) * mp_binom(abs(m), lx - 2 * k)
        coeff += i_fact * k_sum
    coeff *= mp_sqrt(
        mp_fact(2 * lx)
        * mp_fact(2 * ly)
        * mp_fact(2 * lz)
        * lfact
        * mp_fact(lmm)
        / (
            mp_fact(2 * l)
            * mp_fact(lx)
            * mp_fact(ly)
            * mp_fact(lz)
            * mp_fact(l + abs(m))
        )
    )
    coeff *= 1 / (2**l * lfact)

    # print(f"{lx=}, {ly=}, {lz=}, {m=}, {coeff}")
    return coeff


@functools.cache
def cart2sph_coeffs_for(
    l: int,
    real: bool = True,
    zero_small: bool = True,
    zero_thresh: float = ZERO_THRESH,
    use_mp=False,
) -> np.ndarray:
    if use_mp:
        cart2sph_coeff_func = mp_cart2sph_coeff
        sqrt_func = mp_sqrt
        # dtype = np.longcomplex
    else:
        cart2sph_coeff_func = cart2sph_coeff
        sqrt_func = sqrt
        # dtype = complex
    all_coeffs = list()
    for lx, ly, lz in canonical_order(l):
        coeffs = list()
        for m in range(-l, l + 1):
            coeff = cart2sph_coeff_func(lx, ly, lz, m)
            # Form linear combination for m != 0 to get real spherical orbitals.
            if real and m != 0:
                coeff_minus = cart2sph_coeff_func(lx, ly, lz, -m)
                sign = 1 if m < 0 else -1
                coeff = (coeff + sign * coeff_minus) / sqrt_func(sign * 2)
            coeffs.append(coeff)
        all_coeffs.append(coeffs)
    C = np.array(all_coeffs, dtype=complex)  # or use dtype=dtype
    if real:
        assert real == ~np.iscomplex(C).all()  # C should be real
        C = np.real(C)
    if zero_small:
        mask = np.abs(C) <= zero_thresh
        C[mask] = 0.0
    # Final shape will be (spherical, Cartesian)
    # Return read-only view; otherwise it would be too easy to mess up the
    # cached result by inplace operations on C.
    C = C.T.view()
    C.flags.writeable = False
    return C


@functools.cache
def expand_sph_quantum_numbers(
    Lm, zero_thresh=ZERO_THRESH, with_lmn_factors=False
) -> Tuple[np.ndarray, List[Tuple[int, int, int]]]:
    """Factors and Cart. angular momentum vectors for given sph. quantum numbers."""
    L, m = Lm
    m += L  # Map m from [-L, L] onto [0, 2*L]
    C = cart2sph_coeffs_for(L)  # shape (spherical, Cartesian)
    # Include lmn-factor that appears when contracted functions are normalized
    # according to 'normalization.norm_cgto_lmn(L).
    if with_lmn_factors:
        lmn_factors = get_lmn_factors(L)
        C = C * lmn_factors[None, :]
    indices, *_ = np.where(np.abs(C[m]) > zero_thresh)
    factors = C[m, indices].copy()

    cart_lmn = [lmn for i, lmn in enumerate(canonical_order(L)) if i in indices]
    return factors, cart_lmn


def cart2sph_coeffs(
    l_max: int,
    **kwargs,
) -> Dict[int, np.ndarray]:
    coeffs = {l: cart2sph_coeffs_for(l, **kwargs) for l in range(l_max + 1)}
    return coeffs


def cart2sph_nlms(l_max: int) -> Dict[int, Tuple[Tuple[int, int, int]]]:
    nlms = dict()
    for l in range(l_max + 1):
        n = l
        nlms[l] = tuple([(n, l, m) for m in range(-l, l + 1)])
    return nlms


def print_coeffs(l: int, C: np.ndarray):
    print("".join([f"{m: >13d}" for m in range(-l, l + 1)]))
    for (lx, ly, lz), line in zip(canonical_order(l), C):
        fmt = "{:>12.4f}" * len(line)
        verb = "".join(["x"] * lx + ["y"] * ly + ["z"] * lz)
        print(f"{verb}: {fmt.format(*line)}")


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("l_max", type=int, default=3)

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    l_max = args.l_max

    Cs = cart2sph_coeffs(l_max)
    for l, C in Cs.items():
        print(f"{l=}")
        print_coeffs(l, C.T)
        print()


if __name__ == "__main__":
    run()
