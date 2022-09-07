#!/usr/bin/env python

# [1] https://doi.org/10.1002/qua.560540202
#     Transformation between Cartesian and pure spherical harmonic Gaussians
#     Schlegel, Frisch, 1995

import argparse
from cmath import sqrt
from math import floor
import sys
from typing import Dict

import numpy as np
from numpy.typing import NDArray
from scipy.special import factorial as fact
from scipy.special import binom

from pysisyphus.wavefunction.helpers import canonical_order


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
    coeff *= 1 / (2 ** l * lfact)
    return coeff


def cart2sph_coeffs_for(l: int, real: bool = True) -> NDArray:
    all_coeffs = list()
    for lx, ly, lz in canonical_order(l):
        coeffs = list()
        for m in range(-l, l + 1):
            coeff = cart2sph_coeff(lx, ly, lz, m)
            # Form linear combination for m != 0 to get real spherical orbitals.
            if real and m != 0:
                coeff_minus = cart2sph_coeff(lx, ly, lz, -m)
                sign = 1 if m < 0 else -1
                coeff = (coeff + sign * coeff_minus) / sqrt(sign * 2)
            coeffs.append(coeff)
        all_coeffs.append(coeffs)
    C = np.array(all_coeffs, dtype=complex)
    if real:
        assert real == ~np.iscomplex(C).all()  # C should be real
        C = np.real(C)
    return C.T  # Return w/ shape (sph., cart.)


def cart2sph_coeffs(l_max: int, **kwargs) -> Dict[int, NDArray]:
    coeffs = {l: cart2sph_coeffs_for(l, **kwargs) for l in range(l_max + 1)}
    return coeffs


def print_coeffs(l: int, C: NDArray):
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
