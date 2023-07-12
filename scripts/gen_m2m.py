#!/usr/bin/env python3

"""
Generate python code to shift multipoles.
See pysisyphus.wavefunction.dma for more information.
"""

import black
import itertools as it

import numpy as np
from scipy.special import binom
import sympy as sym

from pysisyphus.wavefunction.helpers import lm_iter
from pysisyphus.wavefunction.dma import get_index_maps

from sympleints.sym_solid_harmonics import Rlm as Rlm_sym


PREC = 24
PRINT_PREC = 15
_LE_MAX = 2  # Up to quadrupoles
# _LE_MAX = 3  # Up to Octopoles
_TWO_LE_MAX = 2 * _LE_MAX
_SQRT2 = 1 / sym.sqrt(2.0).evalf(PREC)
# Map between 2d (l, m) and 1d-indices.
_INDEX_MAP = get_index_maps(_LE_MAX)


def get_C(lmax: int = _TWO_LE_MAX) -> dict[(int, int), complex]:
    """Factors for complex->real transformation of multipoles.

    Defined below A9 in [3].
    """
    C = dict()
    for ma in range(-lmax, lmax + 1):
        for mb in range(-lmax, lmax + 1):
            key = (ma, mb)
            if key == (0, 0):
                value = 1.0
            elif abs(ma) != abs(mb):
                value = 0.0
            elif ma > 0 and mb > 0 and ma == mb:
                value = (-1) ** ma * _SQRT2
            elif ma > 0 and (-ma == mb):
                value = _SQRT2
            elif ma < 0 and (-ma == mb):
                value = -1j * (-1) ** mb * _SQRT2
            elif ma < 0 and mb < 0 and ma == mb:
                value = 1j * _SQRT2
            else:
                raise Exception("Safeguard")
            C[key] = value
    return C


C = get_C()
# Complex conjugate values of C
CCONJ = {k: np.conj(v) for k, v in C.items()}


def get_W_sym(L_max: int = _LE_MAX):
    x, y, z, r = sym.symbols("x y z r", real=True)
    r2sub = x**2 + y**2 + z**2
    W = dict()
    keys = list()
    for (l, m), (ldash, mdash) in it.product(lm_iter(L_max), lm_iter(L_max)):
        lm_ind = _INDEX_MAP[(l, m)]
        lmdash_ind = _INDEX_MAP[(ldash, mdash)]
        key = (lmdash_ind, lm_ind)
        keys.append(key)
        for m1 in range(-l, l + 1):
            Cmm1 = C[(m, m1)]
            if Cmm1 == 0.0:
                continue
            for m2 in range(-ldash, ldash + 1):
                Cmdashm2 = CCONJ[(mdash, m2)]
                if Cmdashm2 == 0.0:
                    continue
                binom1 = binom(l + m1, ldash + m2)
                binom2 = binom(l - m1, ldash - m2)
                prefact = sym.sqrt(binom1 * binom2).evalf(PREC)
                if prefact.evalf() == 0.0:
                    continue
                for m3 in range(ldash - l, l - ldash + 1):
                    Cprod = CCONJ[(m3, m1 - m2)]
                    if Cprod == 0.0:
                        continue
                    Cprod *= Cmm1 * Cmdashm2
                    # Index order is consistent with eq. (A11) in [3] (W_{l'm', lm}).
                    rlm = Rlm_sym(l - ldash, m3, x, y, z).subs(r**2, r2sub)
                    summand = prefact * Cprod * rlm
                    W.setdefault(key, 0.0)
                    W[key] += summand
    return W, keys


L_max = _LE_MAX
num = (2 * np.arange(L_max + 1) + 1).sum()
# Nonzero expressions and all keys
Wsym, all_keys = get_W_sym(L_max)

# Carry out all possible evaluations
for key, expr in Wsym.items():
    expr = sym.simplify(expr)
    expr = expr.evalf(PREC)
    Wsym[key] = expr

# All non-zero keys
present_keys = set(Wsym.keys())
# Search for keys in all_keys for values that aren't 0.0.
nonzero_keys = [key for key in all_keys if key in present_keys]
# We can select nonzero_keys to get the expressions in the correct order
W_exprs = [Wsym[key] for key in nonzero_keys]

# Common subexpression elimination
repls, reduced = sym.cse(W_exprs, order="none", optimizations="basic")

Wmat = sym.ImmutableSparseMatrix(num, num, Wsym)
# Dummy multipole matrix
Q = sym.MatrixSymbol("Q", 1, num)
# Carry out matrix multiplication of multipoles with M2M operator
res = sym.Matrix(Q @ Wmat)


def expr2str(expr, callback=None):
    expr = sym.nfloat(expr, n=PRINT_PREC)
    expr = expr.subs(1.0, 1)
    line = str(expr)
    # I'm sorry mom! I just wanted to get rid of these multiplications with 1.0!
    line = line.replace("1.0*", "")
    if callback is not None:
        line = callback(line)
    return line


def replQ(line):
    return line.replace("Q[0, ", "Q[")


repls, reduced = sym.cse(list(res), order="none", optimizations="basic")

# Transform replacement expressions to lines
repl_lines = [f"{lhs} = {expr2str(rhs, callback=replQ)}" for lhs, rhs in repls]
# Transform reduced expressions to lines
lines = [expr2str(expr, callback=replQ) for expr in reduced]
# Join & render them
repl_joined = "\n".join(repl_lines)
joined = ",\n".join(lines)
tpl = f"""# {L_max=}, moves {num} spherical multipoles

{repl_joined}

np.array(({joined}))"""
# Format with black
rendered = black.format_str(tpl, mode=black.FileMode(line_length=90))
print()
print(rendered)
