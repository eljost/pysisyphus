#!/usr/bin/env python3

# [1] https://doi.org/10.1063/1.450106
#     Efficient recursive computation of molecular integrals over Cartesian
#     Gaussian functions
#     Obara, Saika, 1986
# [2] https://doi.org/10.1002/9781119019572
#     Molecular Electronic-Structure Theory
#     Helgaker, Jørgensen, Olsen
# [3] https://doi.org/10.1021/acs.jctc.7b00788
#     LIBRETA: Computerized Optimization and Code Synthesis for
#     Electron Repulsion Integral Evaluation
#     Jun Zhang
# [4] https://doi.org/10.1039/B413539C
#     Efficient evaluation of three-center two-electron integrals over Gaussian functions
#     Ahlrichs, 2004
# [5] EVALUATING MANY-ELECTRON MOLECULAR INTEGRALS FOR QUANTUM CHEMISTRY
#     James Christopher Womack, PhD Thesis
# [6] https://doi.org/10.1063/1.4983393
#     Efficient evaluation of three-center Coulomb integrals
#     Samu, Kállay, 2017
# [7] https://arxiv.org/pdf/2210.03192.pdf
#     Memory-Efficient Recursive Evaluation of 3-Center Gaussian Integrals
#     Asadchev, Valeev, 2022


import argparse
from datetime import datetime
import functools
import itertools as it
import os
from pathlib import Path
import random
import string
import sys
import textwrap
import time

from jinja2 import Template
import numpy as np
from sympy import (
    Array,
    cse,
    exp,
    Expr,
    flatten,
    Function,
    IndexedBase,
    Matrix,
    permutedims,
    pi,
    sqrt,
    Symbol,
    symbols,
    tensorcontraction as tc,
    tensorproduct as tp,
)


from sympy.codegen.ast import Assignment
from sympy.printing.numpy import NumPyPrinter
from sympy.printing.c import C99CodePrinter

# from pysisyphus.wavefunction.cart2sph import cart2sph_coeffs

try:
    from pysisyphus.config import L_MAX, L_AUX_MAX
except ModuleNotFoundError:
    L_MAX = 4
    L_AUX_MAX = 5


L_MAP = {
    0: "s",
    1: "p",
    2: "d",
    3: "f",
    4: "g",
    5: "h",
}

KEYS = (
    "cgto",
    "ovlp",
    "dpm",
    "dqpm",
    "qpm",
    "kin",
    "coul",
    "3c2e",
    "3c2e_sph",
)
ONE_THRESH = 1e-14


def make_py_func(repls, reduced, args=None, name=None, doc_str=""):
    if args is None:
        args = list()
    # Generate random name, if no name was supplied
    if name is None:
        name = "func_" + "".join(
            [random.choice(string.ascii_letters) for i in range(8)]
        )

    # This allows using the 'boys' function without producing an error
    print_settings = {
        "allow_unknown_functions": True,
    }
    print_func = NumPyPrinter(print_settings).doprint
    assignments = [Assignment(lhs, rhs) for lhs, rhs in repls]
    py_lines = [print_func(as_) for as_ in assignments]
    return_val = print_func(reduced)

    tpl = Template(
        """
    def {{ name }}({{ args }}):
        {% if doc_str %}
        \"\"\"{{ doc_str }}\"\"\"
        {% endif %}

        {% for line in py_lines %}
        {{ line }}
        {% endfor %}

        # {{ n_return_vals }} item(s)
        return numpy.array({{ return_val }})
    """,
        trim_blocks=True,
        lstrip_blocks=True,
    )

    rendered = textwrap.dedent(
        tpl.render(
            name=name,
            args=args,
            py_lines=py_lines,
            return_val=return_val,
            n_return_vals=len(reduced),
            doc_str=doc_str,
        )
    ).strip()
    return rendered


def make_c_func(repls, reduced, args=None, name=None, doc_str=""):
    if args is None:
        args = list()
    # Generate random name, if no name was supplied
    if name is None:
        name = "func_" + "".join(
            [random.choice(string.ascii_letters) for i in range(8)]
        )
    arg_strs = list()
    for arg in args:
        if arg.islower():
            arg_str = f"double {arg}"
        elif arg.isupper():
            arg_str = f"double {arg}[3]"
        else:
            raise Exception
        arg_strs.append(arg_str)
    args_str = ", ".join(arg_strs)

    # This allows using the 'boys' function without producing an error
    print_settings = {
        "allow_unknown_functions": True,
        # Without disabling contract some expressions will raise ValueError.
        "contract": False,
    }
    print_func = C99CodePrinter(print_settings).doprint
    assignments = [Assignment(lhs, rhs) for lhs, rhs in repls]
    repl_lines = [print_func(as_) for as_ in assignments]
    res_lines = [print_func(red) for red in reduced]
    res_len = len(reduced)

    signature = f"double * {name}({args_str})"

    tpl = Template(
        """
    {{ signature }} {
        {% if doc_str %}
        /* {{ doc_str }} */
        {% endif %}

        static double {{ res_name }}[{{ res_len}}];

        {% for line in repl_lines %}
        const double {{ line }}
        {% endfor %}

        {% for rhs in res_lines %}
        {{ res_name }}[{{ loop.index0}}] = {{ rhs }};
        {% endfor %}

        return {{ res_name }};
    }
    """,
        trim_blocks=True,
        lstrip_blocks=True,
    )

    rendered = textwrap.dedent(
        tpl.render(
            signature=signature,
            res_name="results",
            res_len=res_len,
            repl_lines=repl_lines,
            res_lines=res_lines,  # c_lines=c_lines,
            reduced=reduced,
            doc_str=doc_str,
            args_str=args_str,
        )
    ).strip()
    return rendered, signature


def canonical_order(L):
    inds = list()
    for i in range(L + 1):
        l = L - i
        for n in range(i + 1):
            m = i - n
            inds.append((l, m, n))
    return inds


def shell_iter(Ls):
    """Iterator over cartesian product of L values in Ls."""
    return it.product(*[canonical_order(L) for L in Ls])


class CartGTO3d(Function):
    """3D Cartesian Gaussian function; not normalized."""

    @classmethod
    @functools.cache
    def eval(cls, i, j, k, a, Xa, Ya, Za):
        Xa2 = Xa**2
        Ya2 = Ya**2
        Za2 = Za**2
        return (Xa**i) * (Ya**j) * (Za**k) * exp(-a * (Xa2 + Ya2 + Za2))


class CartGTOShell(Function):
    @classmethod
    def eval(cls, La_tot, a, Xa, Ya, Za):
        exprs = [CartGTO3d(*La, a, Xa, Ya, Za) for La, in shell_iter((La_tot,))]
        # print(CartGTO3d.eval.cache_info())
        return exprs


class TwoCenter1d(Expr):
    def __new__(self, a, A, b, B, C=None):
        self.a = a
        self.A = A
        self.b = b
        self.B = B
        self.C = C
        return super().__new__(self, a, A, b, B, C)

    @property
    def p(self):
        """Total exponent p."""
        return self.a + self.b

    @property
    def mu(self):
        """Reduced exponent mu."""
        return self.a * self.b / self.p

    @property
    def AB(self):
        """Relative coordinate/Gaussian separation X_AB."""
        return self.A - self.B

    @property
    def P(self):
        """Center-of-charge coordinate."""
        return (self.a * self.A + self.b * self.B) / self.p

    @property
    def PA(self):
        """Relative coordinate/Gaussian separation X_PA."""
        return self.P - self.A

    @property
    def PB(self):
        """Relative coordinate/Gaussian separation X_PB."""
        return self.P - self.B

    @property
    def PC(self):
        """Relative coordinate/Gaussian separation X_PC."""
        return self.P - self.C

    @property
    def K(self):
        return exp(-self.mu * self.AB**2)


class Multipole1d(TwoCenter1d):
    """1d multipole-moment integral of order 'e', between primitive 1d Gaussians
    Ga = G_i(a, r, A) and Gb = G_j(b, r, B) with Cartesian quantum number i and j,
    exponents a and b, centered at A (B). The origin of the multipole expansion is
    at C.
    """

    @functools.cache
    def eval(self, i, j, e):
        ang_moms = (i, j, e)
        if any([_ < 0 for _ in ang_moms]):
            return 0

        recur = self.eval

        def vrr(i, j, e, X):
            return X * recur(i, j, e) + 1 / (2 * self.p) * (
                i * recur(i - 1, j, e) + j * recur(i, j - 1, e) + e * recur(i, j, e - 1)
            )

        # Base case
        if all([_ == 0 for _ in ang_moms]):
            return sqrt(pi / self.p) * self.K
        # Decrement i
        elif i.is_positive:
            return vrr(i - 1, j, e, self.PA)
        # Decrement j
        elif j.is_positive:
            return vrr(i, j - 1, e, self.PB)
        elif e.is_positive:
            return vrr(i, j, e - 1, self.PC)


class Multipole3d(Function):
    @classmethod
    def eval(cls, La, Lb, a, b, A, B, Le, C):
        x, y, z = [
            Multipole1d(a, A[i], b, B[i], C[i]).eval(La[i], Lb[i], Le[i])
            for i in range(3)
        ]
        return x * y * z


class Multipole3dShell(Function):
    @classmethod
    def eval(cls, La_tot, Lb_tot, a, b, A, B, Le_tot=0, C=(0.0, 0.0, 0.0)):
        exprs = [
            Multipole3d(La, Lb, a, b, A, B, Le, C)
            for Le, La, Lb in shell_iter((Le_tot, La_tot, Lb_tot))
        ]
        # print(Multipole1d.eval.cache_info())
        return exprs


class DiagQuadrupole3dShell(Function):
    @classmethod
    def eval(cls, La_tot, Lb_tot, a, b, A, B, C=(0.0, 0.0, 0.0)):
        exprs = list()
        for Le in ((2, 0, 0), (0, 2, 0), (0, 0, 2)):
            for La, Lb in shell_iter((La_tot, Lb_tot)):
                exprs.append(Multipole3d(La, Lb, a, b, A, B, Le, C))
        # print(DiagQuadrupole3dShell.eval.cache_info())
        return exprs


class Overlap1d(TwoCenter1d):
    def eval(self, i, j):
        return Multipole1d(self.a, self.A, self.b, self.B).eval(i, j, 0)


class Overlap3dShell(Function):
    """Whole shell of 3d overlap integrals for a given L pair
    over primitive gaussians."""

    @classmethod
    def eval(cls, La_tot, Lb_tot, a, b, A, B):
        exprs = Multipole3dShell(La_tot, Lb_tot, a, b, A, B)
        # print(Multipole1d.eval.cache_info())
        return exprs


class Kinetic1d(TwoCenter1d):
    @functools.cache
    def eval(self, i, j):
        if i < 0 or j < 0:
            return 0

        recur = self.eval

        def recur_rel(i, j, X):
            return X * recur(i, j) + 1 / (2 * self.p) * (
                i * recur(i - 1, j) + j * recur(i, j - 1)
            )

        def recur_ovlp(i, j):
            return Overlap1d(self.a, self.A, self.b, self.B).eval(i, j)

        # Base case
        if i == 0 and j == 0:
            return (
                self.a - 2 * self.a**2 * (self.PA**2 + 1 / (2 * self.p))
            ) * recur_ovlp(i, j)
        # Decrement i
        elif i > 0:
            # Eq. (9.3.41)
            return recur_rel(i - 1, j, self.PA) + self.b / self.p * (
                2 * self.a * recur_ovlp(i, j) - i * recur_ovlp(i - 2, j)
            )
        # Decrement j
        elif j > 0:
            # Eq. (9.3.41)
            return recur_rel(i, j - 1, self.PB) + self.a / self.p * (
                2 * self.b * recur_ovlp(i, j) - j * recur_ovlp(i, j - 2)
            )


class Kinetic3d(Function):
    @classmethod
    def eval(cls, La, Lb, a, b, A, B):
        Tx, Ty, Tz = [Kinetic1d(a, A[i], b, B[i]).eval(La[i], Lb[i]) for i in range(3)]
        Sx, Sy, Sz = [Overlap1d(a, A[i], b, B[i]).eval(La[i], Lb[i]) for i in range(3)]
        return Tx * Sy * Sz + Sx * Ty * Sz + Sx * Sy * Tz


class Kinetic3dShell(Function):
    @classmethod
    def eval(cls, La_tot, Lb_tot, a, b, A, B):
        exprs = [
            Kinetic3d(La, Lb, a, b, A, B) for La, Lb in shell_iter((La_tot, Lb_tot))
        ]
        return exprs


# Placeholder for the Boys-function. The actual Boys-function will be
# imported in the generated python module.
boys = Function("boys")


class Coulomb(TwoCenter1d):
    """Nucleus at C."""

    @functools.cache
    def eval(self, i, k, m, j, l, n, N):
        ang_moms = (i, k, m, j, l, n)
        if any([am < 0 for am in ang_moms]):
            return 0

        def recur(N, *inds):
            """Simple wrapper to pass all required arguments."""
            return self.eval(*inds, N)  # , a, b, A, B, C)

        def decr(to_decr, decr_ind):
            one = np.zeros(3, dtype=int)
            one[decr_ind] = 1

            def Ni(inds):
                return inds[decr_ind]

            bra = np.array((i, k, m), dtype=int)
            ket = np.array((j, l, n), dtype=int)

            if to_decr == "bra":
                X = self.PA
                bra[decr_ind] -= 1
            else:
                X = self.PB
                ket[decr_ind] -= 1

            bra_decr = bra - one
            ket_decr = ket - one

            return (
                X[decr_ind] * recur(N, *bra, *ket)
                - self.PC[decr_ind] * recur(N + 1, *bra, *ket)
                + 1
                / (2 * self.p)
                * Ni(bra)
                * (recur(N, *bra_decr, *ket) - recur(N + 1, *bra_decr, *ket))
                + 1
                / (2 * self.p)
                * Ni(ket)
                * (recur(N, *bra, *ket_decr) - recur(N + 1, *bra, *ket_decr))
            )

        # Base case
        if all([am == 0 for am in ang_moms]):
            K = exp(-self.mu * self.AB.dot(self.AB))
            return 2 * pi / self.p * K * boys(N, self.p * self.PC.dot(self.PC))
        elif i > 0:
            return decr("bra", 0)
        elif j > 0:
            return decr("ket", 0)
        elif k > 0:
            return decr("bra", 1)
        elif l > 0:
            return decr("ket", 1)
        elif m > 0:
            return decr("bra", 2)
        elif n > 0:
            return decr("ket", 2)


class CoulombShell(Function):
    @classmethod
    def eval(cls, La_tot, Lb_tot, a, b, A, B, C=(0.0, 0.0, 0.0)):
        exprs = [
            # Coulomb(*La, *Lb, 0, a, b, A, B, C)
            Coulomb(a, A, b, B, C).eval(*La, *Lb, 0)
            for La, Lb in shell_iter((La_tot, Lb_tot))
        ]
        # print(Coulomb.eval.cache_info())
        return exprs


class TwoCenterTwoElectron(Function):
    @classmethod
    @functools.cache
    def eval(cls, ia, ja, ka, ib, jb, kb, N, a, b, A, B):
        ang_moms = np.array((ia, ja, ka, ib, jb, kb), dtype=int)
        ang_moms2d = ang_moms.reshape(-1, 3)
        if any([am < 0 for am in ang_moms]):
            return 0

        p = a + b
        P = (a * A + b * B) / p
        mu = (a * b) / p

        def recur(N, *inds):
            return cls(*inds, N, a, b, A, B)

        def vrr(bra_or_ket, cart_ind):
            assert bra_or_ket in ("bra", "ket")

            if bra_or_ket == "bra":
                ind1 = 0
                ind2 = 1
            else:
                ind1 = 1
                ind2 = 0
            exps = (a, b)
            X = (P - A, P - B)[ind1][cart_ind]

            l1 = ang_moms2d[ind1, cart_ind] - 1
            exp1 = exps[ind1]
            decr1 = ang_moms2d.copy()
            decr1[ind1, cart_ind] -= 1
            decr11 = decr1.copy()
            decr11[ind1, cart_ind] -= 1

            l2 = ang_moms2d[ind2, cart_ind]
            exp2 = exps[ind2]
            decr12 = ang_moms2d.copy()
            decr12[ind1, cart_ind] -= 1
            decr12[ind2, cart_ind] -= 1

            decr1 = decr1.flatten()
            decr11 = decr11.flatten()
            decr12 = decr12.flatten()

            return (
                X * recur(N + 1, *decr1)
                + l1
                / (2 * exp1)
                * (recur(N, *decr11) - exp2 / p * recur(N + 1, *decr11))
                + l2 / (2 * p) * recur(N + 1, *decr12)
            )

        # vrr_bra = functools.partial(vrr, bra_or_ket="bra")
        vrr_bra = functools.partial(vrr, "bra")
        vrr_ket = functools.partial(vrr, "ket")

        # Base case
        if (ang_moms == 0).all():
            AB = A - B
            return 2 * pi**2.5 / sqrt(p) / (a * b) * boys(N, mu * AB.dot(AB))
        elif ia > 0:
            return vrr_bra(0)
        elif ja > 0:
            return vrr_bra(1)
        elif ka > 0:
            return vrr_bra(2)
        elif ib > 0:
            return vrr_ket(0)
        elif jb > 0:
            return vrr_ket(1)
        elif kb > 0:
            return vrr_ket(2)


class TwoCenterTwoElectronShell(Function):
    @classmethod
    def eval(cls, La_tot, Lb_tot, a, b, A, B):
        exprs = [
            TwoCenterTwoElectron(*La, *Lb, 0, a, b, A, B)
            for La, Lb in shell_iter((La_tot, Lb_tot))
        ]
        # print(TwoCenterTwoElectron.eval.cache_info())
        return exprs


class ThreeCenterTwoElectronBase(Function):
    """
    https://pubs.rsc.org/en/content/articlelanding/2004/CP/b413539c

    There is an error in the base case (00|0). One must divide by
    sqrt(eta + gamma), not multiply.
    """

    @classmethod
    @functools.cache
    def eval(cls, ia, ja, ka, ib, jb, kb, ic, jc, kc, N, a, b, c, A, B, C):
        ang_moms = np.array((ia, ja, ka, ib, jb, kb, ic, jc, kc), dtype=int)
        ang_moms2d = ang_moms.reshape(-1, 3)
        if any([am < 0 for am in ang_moms]):
            return 0

        p = a + b
        P = (a * A + b * B) / p
        mu = (a * b) / p
        X_PC = P - C

        rho = p * c / (p + c)

        def recur(N, *inds):
            """Simple wrapper to pass all required arguments.

            Here we don't use ThreeCenterTwoElectronBase, but the derived classes,
            that provide the 'aux_vrr' attribute, to distinguish between the different
            vertical recursion relations."""
            return cls(*inds, N, a, b, c, A, B, C)

        def recur_hrr(cart_ind):
            """Horizontal recursion relation to transfer angular momentum in bra.

            (a, b+1_i|c) = (a + 1_i,b|c) + X_AB (ab|c)
            """
            assert N == 0
            incr_ang_moms = ang_moms2d.copy()
            incr_ang_moms[0, cart_ind] += 1  # Increase in left orbital
            incr_ang_moms[1, cart_ind] -= 1  # Decrease in right orbital

            decr_ang_moms = ang_moms2d.copy()
            decr_ang_moms[1, cart_ind] -= 1  # Decrease in right orbital

            incr_ang_moms = incr_ang_moms.flatten()
            decr_ang_moms = decr_ang_moms.flatten()

            AB_dir = (A - B)[cart_ind]

            return recur(N, *incr_ang_moms) + AB_dir * recur(N, *decr_ang_moms)

        def recur_vrr(cart_ind):
            assert (ib, jb, kb) == (0, 0, 0)
            assert (ic, jc, kc) == (0, 0, 0)

            decr_a = ang_moms2d.copy()
            decr_a[0, cart_ind] -= 1  # Decrease in bra left orbital
            decr_aa = decr_a.copy()
            decr_aa[0, cart_ind] -= 1  # Decrease in bra left orbital again

            PA_dir = (P - A)[cart_ind]
            PC_dir = (P - C)[cart_ind]
            ai = (ia, ja, ka)[cart_ind] - 1
            _2p = 2 * p

            decr_a = decr_a.flatten()
            decr_aa = decr_aa.flatten()

            return (
                PA_dir * recur(N, *decr_a)
                - rho / p * PC_dir * recur(N + 1, *decr_a)
                + ai / _2p * (recur(N, *decr_aa) - rho / p * recur(N + 1, *decr_aa))
            )

        def recur_vrr_aux(cart_ind):
            decr_c = ang_moms2d.copy()
            decr_c[2, cart_ind] -= 1

            decr_cc = decr_c.copy()
            decr_cc[2, cart_ind] -= 1

            decr_ac = decr_c.copy()
            decr_ac[0, cart_ind] -= 1

            decr_bc = decr_c.copy()
            decr_bc[1, cart_ind] -= 1

            PC_dir = (P - C)[cart_ind]
            la = (ia, ja, ka)[cart_ind]
            lb = (ib, jb, kb)[cart_ind]
            lc = (ic, jc, kc)[cart_ind] - 1

            decr_c = decr_c.flatten()
            decr_cc = decr_cc.flatten()
            decr_ac = decr_ac.flatten()
            decr_bc = decr_bc.flatten()

            return (
                p / (p + c) * PC_dir * recur(N + 1, *decr_c)
                + lc
                / (2 * c)
                * (recur(N, *decr_cc) - p / (p + c) * recur(N + 1, *decr_cc))
                + la / (2 * (p + c)) * recur(N + 1, *decr_ac)
                + lb / (2 * (p + c)) * recur(N + 1, *decr_bc)
            )

        def recur_vrr_aux_sph(cart_ind):
            assert (ib, jb, kb) == (0, 0, 0)
            decr_c = ang_moms2d.copy()
            decr_c[2, cart_ind] -= 1
            decr_ac = decr_c.copy()
            decr_ac[0, cart_ind] -= 1
            La = (ia, ja, ka)[cart_ind]
            PC_dir = (P - C)[cart_ind]
            return (
                rho
                / c
                * (
                    PC_dir * recur(N + 1, *decr_c.flatten())
                    + La / (2 * p) * recur(N + 1, *decr_ac.flatten())
                )
            )

        recur_vrr_aux_funcs = {
            "cart": recur_vrr_aux,
            "sph": recur_vrr_aux_sph,
        }
        recur_vrr_aux_func = recur_vrr_aux_funcs[cls.aux_vrr]

        # Base case
        if (ang_moms == 0).all():
            X_AB = A - B
            r2_PC = X_PC.dot(X_PC)
            r2_AB = X_AB.dot(X_AB)
            chi = rho * r2_PC
            K = exp(-mu * r2_AB)
            return 2 * pi**2.5 / sqrt(p + c) / (p * c) * K * boys(N, chi)
        elif ib > 0:
            return recur_hrr(0)
        elif jb > 0:
            return recur_hrr(1)
        elif kb > 0:
            return recur_hrr(2)
        elif ic > 0:
            return recur_vrr_aux_func(0)
        elif jc > 0:
            return recur_vrr_aux_func(1)
        elif kc > 0:
            return recur_vrr_aux_func(2)
        elif ia > 0:
            return recur_vrr(0)
        elif ja > 0:
            return recur_vrr(1)
        elif ka > 0:
            return recur_vrr(2)


class ThreeCenterTwoElectron(ThreeCenterTwoElectronBase):
    aux_vrr = "cart"


class ThreeCenterTwoElectronSph(ThreeCenterTwoElectronBase):
    aux_vrr = "sph"


class ThreeCenterTwoElectronShell(Function):
    @classmethod
    def eval(cls, La_tot, Lb_tot, Lc_tot, a, b, c, A, B, C):
        exprs = [
            ThreeCenterTwoElectron(*La, *Lb, *Lc, 0, a, b, c, A, B, C)
            for La, Lb, Lc in shell_iter((La_tot, Lb_tot, Lc_tot))
        ]
        # print(ThreeCenterTwoElectron.eval.cache_info())
        return exprs


class ThreeCenterTwoElectronSphShell(Function):
    @classmethod
    def eval(cls, La_tot, Lb_tot, Lc_tot, a, b, c, A, B, C):
        exprs = [
            ThreeCenterTwoElectronSph(*La, *Lb, *Lc, 0, a, b, c, A, B, C)
            for La, Lb, Lc in shell_iter((La_tot, Lb_tot, Lc_tot))
        ]
        # print(ThreeCenterTwoElectron.eval.cache_info())
        return exprs


class ERI(Function):
    """Variables named according to libreta paper [3]."""

    @classmethod
    def eval(
        cls, ia, ja, ka, ib, jb, kb, ic, jc, kc, id_, jd, kd, N, a, b, c, d, A, B, C, D
    ):
        ang_moms = np.array(
            (ia, ja, ka, ib, jb, kb, ic, jc, kc, id_, jd, kd), dtype=int
        )
        ang_moms2d = ang_moms.reshape(-1, 3)

        if any([am < 0 for am in ang_moms]):
            return 0

        AB = A - B
        xi = a + b
        P = (a * A + b * B) / xi  # Eq. (5)

        CD = C - D
        zeta = c + d
        Q = (c * C + d * D) / zeta  # Eq. (6)

        theta = xi * zeta / (xi + zeta)  # Eq. (7)

        def recur(N, *ang_moms):
            """Simple wrapper to pass all required arguments."""
            return ERI(*ang_moms, N, a, b, c, d, A, B, C, D)

        def recur_hrr(bra_or_ket, cart_ind):
            """Horizontal recursion relation to transfer angmoms in bra or ket.

            (a, b+1_x|cd) = (a + 1_x,b|cd) + X_AB (ab|cd)
            (ab|c, d + 1_x) = (ab|c + 1_x, d) + X_CD (ab|cd)
            """
            if bra_or_ket == "bra":
                XYZ = AB
                incr_ind = 0
            else:
                XYZ = CD
                incr_ind = 2

            decr_ind = incr_ind + 1
            incr_ang_moms = ang_moms2d.copy()
            incr_ang_moms[incr_ind, cart_ind] += 1  # Increase in left orbital
            incr_ang_moms[decr_ind, cart_ind] -= 1  # Decrease in right orbital
            incr_ang_moms = incr_ang_moms.flatten()
            decr_ang_moms = ang_moms2d.copy()
            decr_ang_moms[decr_ind, cart_ind] -= 1
            decr_ang_moms = decr_ang_moms.flatten()

            return recur(N, *incr_ang_moms) + XYZ[cart_ind] * recur(N, *decr_ang_moms)

        def recur_vrr(bra_or_ket, cart_ind):
            assert (ib, jb, kb) == (0, 0, 0)
            assert (id_, jd, kd) == (0, 0, 0)

            if bra_or_ket == "bra":
                XYZ = P - A
                decr_ind = 0
                decr_also_ind = 2
                exp1, exp2 = zeta, xi
                am1 = (ia, ja, ka)[cart_ind] - 1
                am2 = (ic, jc, kc)[cart_ind]
                sign = -1
            else:
                XYZ = Q - C
                decr_ind = 2
                decr_also_ind = 0
                exp1, exp2 = xi, zeta
                am1 = (ic, jc, kc)[cart_ind] - 1
                am2 = (ia, ja, ka)[cart_ind]
                sign = 1

            decr_ang_moms = ang_moms2d.copy()
            decr_ang_moms[decr_ind, cart_ind] -= 1
            decr_ang_moms = decr_ang_moms.flatten()
            decr2_ang_moms = ang_moms2d.copy()
            decr2_ang_moms[decr_ind, cart_ind] -= 2  # Further decrease angular momentum
            decr2_ang_moms = decr2_ang_moms.flatten()
            decr_also_ang_moms = ang_moms2d.copy()
            decr_also_ang_moms[decr_also_ind, cart_ind] -= 1
            decr_also_ang_moms = decr_also_ang_moms.flatten()

            denom = xi + zeta
            quot = exp1 / denom
            PQ = P - Q

            return (
                XYZ[cart_ind] * recur(N, *decr_ang_moms)
                + sign * quot * PQ[cart_ind] * recur(N + 1, *decr_ang_moms)
                + am1
                / (2 * exp2)
                * (recur(N, *decr2_ang_moms) - quot * recur(N + 1, *decr2_ang_moms))
                + am2 / (2 * denom) * recur(N + 1, *decr_also_ang_moms)
            )

        # Base case
        if all([am == 0 for am in ang_moms]):
            RAB2 = AB.dot(AB)
            RCD2 = CD.dot(CD)
            Kab = exp(-(a * b) / xi * RAB2)
            Kcd = exp(-(c * d) / zeta * RCD2)

            PQ = P - Q
            RPQ2 = PQ.dot(PQ)

            return (
                2
                * pi ** (5 / 2)
                / (xi * zeta + sqrt(xi + zeta))
                * Kab
                * Kcd
                * boys(N, theta * RPQ2)
            )  # Eq. (8) [ss|ss]^N
        #
        # Horizontal recursion relation to reduce angular momentum.
        #
        # Bra, basis func. with exp. b, centered at B with Lb_tot = (ib + jb + kb)
        elif ib > 0:
            return recur_hrr("bra", 0)
        elif jb > 0:
            return recur_hrr("bra", 1)
        elif kb > 0:
            return recur_hrr("bra", 2)
        # Ket, basis func. with exp. d, centered at D with Ld_tot = (id_ + jd + kd)
        elif id_ > 0:
            return recur_hrr("ket", 0)
        elif jd > 0:
            return recur_hrr("ket", 1)
        elif kd > 0:
            return recur_hrr("ket", 2)
        #
        # Vertical recursion relation to reduce angular momentum.
        #
        # Bra, basis func. with exp. a, centered at A with La_tot = (ia + ja + ka)
        elif ia > 0:
            return recur_vrr("bra", 0)
        elif ja > 0:
            return recur_vrr("bra", 1)
        elif ka > 0:
            return recur_vrr("bra", 2)
        # Ket, basis func. with exp. c, centered at C with Lc_tot = (ic + jc + kc)
        elif ic > 0:
            return recur_vrr("ket", 0)
        elif jc > 0:
            return recur_vrr("ket", 1)
        elif kc > 0:
            return recur_vrr("ket", 1)


class ERIShell(Function):
    @classmethod
    def eval(cls, La_tot, Lb_tot, Lc_tot, Ld_tot, a, b, c, d, A, B, C, D):
        exprs = [
            ERI(*La, *Lb, *Lc, *Ld, 0, a, b, c, d, A, B, C, D)
            for La, Lb, Lc, Ld in shell_iter((La_tot, Lb_tot, Lc_tot, Ld_tot))
        ]
        # print(ERI.eval.cache_info())
        return exprs


class SpinOrbit(Function):
    """1-electron spin-orbit interaction integral for Cartesian direction μ."""

    @classmethod
    @functools.cache
    def eval(cls, i, k, m, j, l, n, N, mu, a, b, A, B, C):
        ang_moms = (i, k, m, j, l, n)
        if any([am < 0 for am in ang_moms]):
            return 0

        p = a + b
        P = (a * A + b * B) / p
        X_PC = P - C

        def coulomb(N, *inds):
            return Coulomb(*inds, N, a, b, A, B, C)

        def recur(N, *inds):
            """Simple wrapper to pass all required arguments."""
            return SpinOrbit(*inds, N, mu, a, b, A, B, C)

        def decr(to_decr, decr_ind):
            one = np.zeros(3, dtype=int)
            one[decr_ind] = 1

            def Ni(inds):
                return inds[decr_ind]

            bra = np.array((i, k, m), dtype=int)
            ket = np.array((j, l, n), dtype=int)

            if to_decr == "bra":
                X = P - A
                bra[decr_ind] -= 1
            else:
                X = P - B
                ket[decr_ind] -= 1

            bra_decr = bra - one
            ket_decr = ket - one

            expr = (
                # The initial part of this recursion is the same as for the 1-electron
                # Coulomb integrals.
                X[decr_ind] * recur(N, *bra, *ket)
                - X_PC[decr_ind] * recur(N + 1, *bra, *ket)
                + 1
                / (2 * p)
                * Ni(bra)
                * (recur(N, *bra_decr, *ket) - recur(N + 1, *bra_decr, *ket))
                + 1
                / (2 * p)
                * Ni(ket)
                * (recur(N, *bra, *ket_decr) - recur(N + 1, *bra, *ket_decr))
                # From here on, the recursion differs.
            )

            if to_decr == "bra":
                X = B - C
                exponent = b
            else:
                X = A - C
                exponent = a
            # Second to line in Eq. (A34) in [1]
            expr += 2 * exponent * Matrix(one).cross(X)[mu] * coulomb(N + 1, *bra, *ket)
            # Last line in Eq. (A34) in [1]
            for k_ in range(
                3
            ):  # Use k_ instead of k, as this is already an angular momentum
                one_k = np.zeros(3, dtype=int)
                one_k[k_] = 1
                one_ind = np.cross(one, one_k)[mu]
                if to_decr == "bra":
                    expr += ket[k_] * one_ind * coulomb(N + 1, *bra, *(ket - one_k))
                else:
                    expr += bra[k_] * one_ind * coulomb(N + 1, *(bra - one_k), *ket)
            return expr

        # Base case
        if all([am == 0 for am in ang_moms]):
            AC = A - C
            BC = B - C
            ACxBC = AC.cross(BC)
            return 4 * p * ACxBC[m] * coulomb(N + 1, 0, 0, 0, 0, 0, 0)
        elif i > 0:
            return decr("bra", 0)
        elif j > 0:
            return decr("ket", 0)
        elif k > 0:
            return decr("bra", 1)
        elif l > 0:
            return decr("ket", 1)
        elif m > 0:
            return decr("bra", 2)
        elif n > 0:
            return decr("ket", 2)


class SpinOrbitShell(Function):
    @classmethod
    def eval(cls, La_tot, Lb_tot, a, b, A, B, C=(0.0, 0.0, 0.0)):
        exprs = [
            SpinOrbit(*La, *Lb, 0, mu, a, b, A, B, C)
            for (La, Lb), mu in it.product(shell_iter((La_tot, Lb_tot)), range(3))
        ]
        # print(SpinOrbit.eval.cache_info())
        return exprs


def get_center(i):
    symbs = [Symbol(str(i) + ind, real=True) for ind in ("x", "y", "z")]
    return Matrix([*symbs]).T  # Return column vector


def get_centers(i):
    symbs = [Symbol(str(i) + ind, real=True) for ind in ("x", "y", "z")]
    return Matrix([*symbs]).T  # Return column vector


def get_map(i, center_i):
    array = IndexedBase(i, shape=3)
    array_map = dict(zip(center_i, array))
    return array, array_map


def cart2spherical(L_tots, exprs):
    assert len(L_tots) > 0

    # Coefficient matrices for Cartesian-to-spherical conversion
    coeffs = [Array(CART2SPH[L]) for L in L_tots]
    cart_shape = [(l + 1) * (l + 2) // 2 for l in L_tots]
    cart = Array(exprs).reshape(*cart_shape)

    sph = tc(tp(coeffs[0], cart), (1, 2))
    if len(L_tots) == 2:
        sph = tc(tp(sph, coeffs[1].transpose()), (1, 2))
    elif len(L_tots) == 3:
        _, Cb, Cc = coeffs
        sph = tc(tp(permutedims(sph, (0, 2, 1)), Cb.transpose()), (2, 3))
        sph = tc(tp(permutedims(sph, (0, 2, 1)), Cc.transpose()), (2, 3))
    else:
        raise Exception(
            "Cartesian -> spherical transformation for 4-center integrals "
            "is not implemented!"
        )

    # Cartesian-to-spherical transformation introduces quite a number of
    # multiplications by 1.0, which are uneccessary. Here, we try to drop
    # some of them by replacing numbers very close to +1.0 with 1.
    sph = sph.replace(lambda n: n.is_Number and (abs(n - 1) <= ONE_THRESH), lambda n: 1)
    # TODO: maybe something along the lines
    # sph = map(lambda expr: expr.evalf(), flatten(sph))
    # is faster?
    return flatten(sph)


def gen_integral_exprs(
    int_func,
    L_maxs,
    kind,
    maps=None,
    sph=False,
):
    if maps is None:
        maps = list()

    ranges = [range(L + 1) for L in L_maxs]

    for L_tots in it.product(*ranges):
        time_str = time.strftime("%H:%M:%S")
        start = datetime.now()
        print(f"{time_str} - Generating {L_tots} {kind}")
        sys.stdout.flush()
        # Generate actual list of expressions.
        exprs = int_func(*L_tots)
        print("\t... generated expressions")
        sys.stdout.flush()
        # if sph:
        # exprs = cart2spherical(L_tots, exprs)
        # print("\t... did Cartesian -> Spherical conversion")
        # sys.stdout.flush()

        # Common subexpression elimination
        repls, reduced = cse(list(exprs), order="none")
        print("\t... did common subexpression elimination")

        # Replacement expressions, used to form the reduced expressions.
        for i, (lhs, rhs) in enumerate(repls):
            rhs = rhs.evalf()
            # Replace occurences of Ax, Ay, Az, ... with A[0], A[1], A[2], ...
            rhs = functools.reduce(lambda rhs, map_: rhs.xreplace(map_), maps, rhs)
            repls[i] = (lhs, rhs)

        # Reduced expression, i.e., the final integrals/expressions.
        for i, red in enumerate(reduced):
            red = red.evalf()
            reduced[i] = functools.reduce(
                lambda red, map_: red.xreplace(map_), maps, red
            )
        # Carry out Cartesian-to-spherical transformation, if requested.
        if sph:
            reduced = cart2spherical(L_tots, reduced)
            print("\t... did Cartesian -> Spherical conversion")
            sys.stdout.flush()

        dur = datetime.now() - start
        print(f"\t... finished in {str(dur)} h")
        sys.stdout.flush()
        yield (repls, reduced), L_tots


def render_py_funcs(exprs_Ls, args, base_name, doc_func, add_imports=None, comment=""):
    if add_imports is None:
        add_imports = ()
    add_imports_str = "\n".join(add_imports)
    if add_imports:
        add_imports_str += "\n\n"

    args = ", ".join((map(str, args)))

    py_funcs = list()
    for (repls, reduced), L_tots in exprs_Ls:
        doc_str = doc_func(L_tots)
        doc_str += "\n\nGenerated code; DO NOT modify by hand!"
        name = base_name + "_" + "".join(str(l) for l in L_tots)
        print(f"Rendering '{name}' ... ", end="")
        start = time.time()
        py_func = make_py_func(repls, reduced, args=args, name=name, doc_str=doc_str)
        dur = time.time() - start
        print(f"finished in {dur: >8.2f} s")
        py_funcs.append(py_func)
    py_funcs_joined = "\n\n".join(py_funcs)

    if comment != "":
        comment = f'"""\n{comment}\n"""\n\n'

    np_tpl = Template(
        "import numpy\n\n{{ add_imports }}{{ comment }}{{ py_funcs }}",
        trim_blocks=True,
        lstrip_blocks=True,
    )
    np_rendered = np_tpl.render(
        comment=comment, add_imports=add_imports_str, py_funcs=py_funcs_joined
    )
    np_rendered = textwrap.dedent(np_rendered)
    try:
        from black import format_str, FileMode

        np_rendered = format_str(np_rendered, mode=FileMode(line_length=90))
    except ModuleNotFoundError:
        print("Skipped formatting with black, as it is not installed!")

    return np_rendered


def render_c_funcs(exprs_Ls, args, base_name, doc_func, add_imports=None, comment=""):
    if add_imports is not None:
        raise Exception("Implement me!")

    arg_strs = [str(arg) for arg in args]

    funcs = list()
    signatures = list()
    for (repls, reduced), L_tots in exprs_Ls:
        doc_str = doc_func(L_tots)
        doc_str += "\n\n\t\tGenerated code; DO NOT modify by hand!"
        doc_str = textwrap.dedent(doc_str)
        name = base_name + "_" + "".join(str(l) for l in L_tots)
        func, signature = make_c_func(
            repls, reduced, args=arg_strs, name=name, doc_str=doc_str
        )
        funcs.append(func)
        signatures.append(signature)
    funcs_joined = "\n\n".join(funcs)

    if comment != "":
        comment = f"/*{comment}*/"

    # Render C files
    c_tpl = Template(
        "#include <math.h>\n\n{{ comment }}\n\n{{ funcs }}",
        trim_blocks=True,
        lstrip_blocks=True,
    )
    c_rendered = c_tpl.render(comment=comment, funcs=funcs_joined)
    c_rendered = textwrap.dedent(c_rendered)
    # Render simple header file
    h_rendered = "\n".join([f"{sig};" for sig in signatures])
    return c_rendered, h_rendered


def write_file(out_dir, name, rendered):
    out_name = out_dir / name
    with open(out_name, "w") as handle:
        handle.write(rendered)
    print(f"Wrote '{out_name}'.")


def write_render(
    ints_Ls,
    args,
    name,
    doc_func,
    out_dir,
    comment="",
    py_kwargs=None,
    c=True,
    c_kwargs=None,
):
    if py_kwargs is None:
        py_kwargs = {}
    if c_kwargs is None:
        c_kwargs = {}
    ints_Ls = list(ints_Ls)
    # Python
    py_rendered = render_py_funcs(
        ints_Ls, args, name, doc_func, comment=comment, **py_kwargs
    )
    write_file(out_dir, f"{name}.py", py_rendered)
    # C
    if c:
        c_rendered, h_rendered = render_c_funcs(
            ints_Ls,
            args,
            name,
            doc_func,
            comment=comment,
            **c_kwargs,
        )
        write_file(out_dir, f"{name}.c", c_rendered)
        write_file(out_dir, f"{name}.h", h_rendered)


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--lmax",
        default=L_MAX,
        type=int,
        help="Generate 1e-integrals up to this maximum angular momentum.",
    )
    parser.add_argument(
        "--lauxmax",
        default=L_AUX_MAX,
        type=int,
        help="Maximum angular moment for integrals using auxiliary functions.",
    )
    parser.add_argument(
        "--write",
        action="store_true",
        help="Write out generated integrals to the current directory, potentially overwriting the present modules.",
    )
    parser.add_argument(
        "--out-dir",
        default="devel_ints",
        help="Directory, where the generated integrals are written.",
    )
    keys_str = f"({', '.join(KEYS)})"
    parser.add_argument(
        "--keys",
        nargs="+",
        help=f"Generate only certain expressions. Possible keys are: {keys_str}. "
        "If not given, all expressions are generated.",
    )

    return parser.parse_args()


def run():
    args = parse_args(sys.argv[1:])

    l_max = args.lmax
    l_aux_max = args.lauxmax
    out_dir = Path(args.out_dir if not args.write else ".")
    keys = args.keys
    if keys is None:
        keys = list()
    try:
        os.mkdir(out_dir)
    except FileExistsError:
        pass

    try:
        global CART2SPH
        CART2SPH = cart2sph_coeffs(max(l_max, l_aux_max), zero_small=True)
    except NameError:
        print("cart2sph_coeffs import is deactivated or pysisyphus is not installed.")

    # Cartesian basis function centers A and B.
    center_A = get_center("A")
    center_B = get_center("B")
    center_C = get_center("C")
    center_D = get_center("D")
    # center_R = get_center("R")
    Xa, Ya, Za = symbols("Xa Ya Za")

    # Orbital exponents a, b, c, d.
    a, b, c, d = symbols("a b c d", real=True)

    # These maps will be used to convert {Ax, Ay, ...} to array quantities
    # in the generated code. This way an iterable/np.ndarray can be used as
    # function argument instead of (Ax, Ay, Az, Bx, By, Bz).
    A, A_map = get_map("A", center_A)
    B, B_map = get_map("B", center_B)
    C, C_map = get_map("C", center_C)
    D, D_map = get_map("D", center_D)

    boys_import = ("from pysisyphus.wavefunction.ints.boys import boys",)

    #################
    # Cartesian GTO #
    #################

    def cart_gto():
        def cart_gto_doc_func(L_tot):
            (La_tot,) = L_tot
            shell_a = L_MAP[La_tot]
            return (
                f"3D Cartesian {shell_a}-Gaussian shell.\n\n"
                "Exponent a, centered at A, evaluated at (Xa, Ya, Za) + A."
            )

        # This code can evaluate multiple points at a time
        cart_gto_Ls = gen_integral_exprs(
            lambda La_tot: CartGTOShell(La_tot, a, Xa, Ya, Za),
            (l_max,),
            "cart_gto",
        )
        cart_gto_rendered = render_py_funcs(
            cart_gto_Ls, (a, Xa, Ya, Za), "cart_gto3d", cart_gto_doc_func
        )
        write_file(out_dir, "gto3d.py", cart_gto_rendered)
        print()

    #####################
    # Overlap integrals #
    #####################

    def overlap():
        def ovlp_doc_func(L_tots):
            La_tot, Lb_tot = L_tots
            shell_a = L_MAP[La_tot]
            shell_b = L_MAP[Lb_tot]
            return f"Cartesian 3D ({shell_a}{shell_b}) overlap integral."

        ovlp_ints_Ls = gen_integral_exprs(
            lambda La_tot, Lb_tot: Overlap3dShell(
                La_tot, Lb_tot, a, b, center_A, center_B
            ),
            (l_max, l_max),
            "overlap",
            (A_map, B_map),
        )
        write_render(
            ovlp_ints_Ls, (a, A, b, B), "ovlp3d", ovlp_doc_func, out_dir, c=True
        )
        print()

    ###########################
    # Dipole moment integrals #
    ###########################

    def dipole():
        def dipole_doc_func(L_tots):
            La_tot, Lb_tot = L_tots
            shell_a = L_MAP[La_tot]
            shell_b = L_MAP[Lb_tot]
            return (
                f"Cartesian 3D ({shell_a}{shell_b}) dipole moment integrals.\n"
                "The origin is at C."
            )

        dipole_comment = """
        Dipole integrals are given in the order:
        for bf_a in basis_functions_a:
            for bf_b in basis_functions_b:
                for cart_dir in (x, y, z):
                    dipole_integrals(bf_a, bf_b, cart_dir)

        So for <s_a|μ|s_b> it will be:

            <s_a|x|s_b>
            <s_a|y|s_b>
            <s_a|z|s_b>
        """

        dipole_ints_Ls = gen_integral_exprs(
            lambda La_tot, Lb_tot: Multipole3dShell(
                La_tot, Lb_tot, a, b, center_A, center_B, 1, center_C
            ),
            (l_max, l_max),
            "dipole moment",
            (A_map, B_map, C_map),
        )
        write_render(
            dipole_ints_Ls,
            (a, A, b, B, C),
            "dipole3d",
            dipole_doc_func,
            out_dir,
            comment=dipole_comment,
            c=True,
        )
        print()

    ###########################################
    # Diagonal of quadrupole moment integrals #
    ###########################################

    def diag_quadrupole():
        def diag_quadrupole_doc_func(L_tots):
            La_tot, Lb_tot = L_tots
            shell_a = L_MAP[La_tot]
            shell_b = L_MAP[Lb_tot]
            return (
                f"Cartesian 3D ({shell_a}{shell_b}) quadrupole moment integrals\n"
                "for operators x², y² and z². The origin is at C."
            )

        diag_quadrupole_comment = """
        Diagonal of the quadrupole moment matrix with operators x², y², z².

        for rr in (xx, yy, zz):
            for bf_a in basis_functions_a:
                for bf_b in basis_functions_b:
                        quadrupole_integrals(bf_a, bf_b, rr)
        """

        diag_quadrupole_ints_Ls = gen_integral_exprs(
            lambda La_tot, Lb_tot: DiagQuadrupole3dShell(
                La_tot, Lb_tot, a, b, center_A, center_B, center_C
            ),
            (l_max, l_max),
            "diag quadrupole moment",
            (A_map, B_map, C_map),
        )
        write_render(
            diag_quadrupole_ints_Ls,
            (a, A, b, B, C),
            "diag_quadrupole3d",
            diag_quadrupole_doc_func,
            out_dir,
            comment=diag_quadrupole_comment,
            c=True,
        )
        print()

    ###############################
    # Quadrupole moment integrals #
    ###############################

    def quadrupole():
        def quadrupole_doc_func(L_tots):
            La_tot, Lb_tot = L_tots
            shell_a = L_MAP[La_tot]
            shell_b = L_MAP[Lb_tot]
            return (
                f"Cartesian 3D ({shell_a}{shell_b}) quadrupole moment integrals.\n"
                "The origin is at C."
            )

        quadrupole_comment = """
        Quadrupole integrals contain the upper triangular part of the symmetric
        3x3 quadrupole matrix.

        / xx xy xz \\
        |    yy yz |
        \       zz /
        """

        quadrupole_ints_Ls = gen_integral_exprs(
            lambda La_tot, Lb_tot: Multipole3dShell(
                La_tot, Lb_tot, a, b, center_A, center_B, 2, center_C
            ),
            (l_max, l_max),
            "quadrupole moment",
            (A_map, B_map, C_map),
        )

        write_render(
            quadrupole_ints_Ls,
            (a, A, b, B, C),
            "quadrupole3d",
            quadrupole_doc_func,
            out_dir,
            comment=quadrupole_comment,
            c=True,
        )
        print()

    ############################
    # Kinetic energy integrals #
    ############################

    def kinetic():
        def kinetic_doc_func(L_tots):
            La_tot, Lb_tot = L_tots
            shell_a = L_MAP[La_tot]
            shell_b = L_MAP[Lb_tot]
            return f"Cartesian 3D ({shell_a}{shell_b}) kinetic energy integral."

        kinetic_ints_Ls = gen_integral_exprs(
            lambda La_tot, Lb_tot: Kinetic3dShell(
                La_tot, Lb_tot, a, b, center_A, center_B
            ),
            (l_max, l_max),
            "kinetic",
            (A_map, B_map),
        )
        write_render(
            kinetic_ints_Ls,
            (a, A, b, B),
            "kinetic3d",
            kinetic_doc_func,
            out_dir,
            c=True,
        )
        print()

    #########################
    # 1el Coulomb Integrals #
    #########################

    def coulomb():
        def coulomb_doc_func(L_tots):
            La_tot, Lb_tot = L_tots
            shell_a = L_MAP[La_tot]
            shell_b = L_MAP[Lb_tot]
            return f"Cartesian ({shell_a}{shell_b}) 1-electron Coulomb integral."

        coulomb_ints_Ls = gen_integral_exprs(
            lambda La_tot, Lb_tot: CoulombShell(
                La_tot, Lb_tot, a, b, center_A, center_B, center_C
            ),
            (l_max, l_max),
            "coulomb3d",
            (A_map, B_map, C_map),
        )
        coulomb_rendered = render_py_funcs(
            coulomb_ints_Ls,
            (a, A, b, B, C),
            "coulomb3d",
            coulomb_doc_func,
            add_imports=boys_import,
        )
        # TODO: handle add_args and Boys function in C
        write_file(out_dir, "coulomb3d.py", coulomb_rendered)
        print()

    ###############################################
    # Two-center two-electron repulsion integrals #
    ###############################################

    def _2center2electron():
        def _2center2el_doc_func(L_tots):
            La_tot, Lb_tot = L_tots
            shell_a = L_MAP[La_tot]
            shell_b = L_MAP[Lb_tot]
            return (
                f"Cartesian ({shell_a}|{shell_b}) "
                "two-center two-electron repulsion integral."
            )

        _2center2el_ints_Ls = gen_integral_exprs(
            lambda La_tot, Lb_tot: TwoCenterTwoElectronShell(
                La_tot,
                Lb_tot,
                a,
                b,
                center_A,
                center_B,
            ),
            (l_aux_max, l_aux_max),
            "_2center2el3d",
            (A_map, B_map),
        )
        write_render(
            _2center2el_ints_Ls,
            (a, A, b, B),
            "_2center2el3d",
            _2center2el_doc_func,
            out_dir,
            c=False,
            py_kwargs={"add_imports": boys_import},
        )
        print()

    #################################################
    # Three-center two-electron repulsion integrals #
    #################################################

    def _3center2electron():
        def _3center2el_doc_func(L_tots):
            La_tot, Lb_tot, Lc_tot = L_tots
            shell_a = L_MAP[La_tot]
            shell_b = L_MAP[Lb_tot]
            shell_c = L_MAP[Lc_tot]
            return (
                f"Cartesian ({shell_a}{shell_b}|{shell_c}) "
                "three-center two-electron repulsion integral."
            )

        _3center2el_ints_Ls = gen_integral_exprs(
            lambda La_tot, Lb_tot, Lc_tot: ThreeCenterTwoElectronShell(
                La_tot, Lb_tot, Lc_tot, a, b, c, center_A, center_B, center_C
            ),
            (l_max, l_max, l_aux_max),
            "_3center2el3d",
            (A_map, B_map, C_map),
        )
        write_render(
            _3center2el_ints_Ls,
            (a, A, b, B, c, C),
            "_3center2el3d",
            _3center2el_doc_func,
            out_dir,
            c=False,
            py_kwargs={"add_imports": boys_import},
        )
        print()

    def _3center2electron_sph():
        def _3center2el_doc_func(L_tots):
            La_tot, Lb_tot, Lc_tot = L_tots
            shell_a = L_MAP[La_tot]
            shell_b = L_MAP[Lb_tot]
            shell_c = L_MAP[Lc_tot]
            return (
                f"Cartesian ({shell_a}{shell_b}|{shell_c}) "
                "three-center two-electron repulsion integral for conversion to "
                "spherical harmonics.\nUses Ahlrichs' vertical recursion relation, "
                "that leaves out some terms, that cancel\nwhen convertig to "
                "spherical harmonics."
            )

        _3center2el_ints_Ls = gen_integral_exprs(
            lambda La_tot, Lb_tot, Lc_tot: ThreeCenterTwoElectronSphShell(
                La_tot, Lb_tot, Lc_tot, a, b, c, center_A, center_B, center_C
            ),
            (l_max, l_max, l_aux_max),
            "_3center2el3d_sph",
            (A_map, B_map, C_map),
            sph=False,
        )
        write_render(
            _3center2el_ints_Ls,
            (a, A, b, B, c, C),
            "_3center2el3d_sph",
            _3center2el_doc_func,
            out_dir,
            c=False,
            py_kwargs={"add_imports": boys_import},
        )
        print()

    ####################################
    # Spin-orbit interaction integrals #
    ####################################

    # def so1el_doc_func(L_tots):
    # La_tot, Lb_tot = L_tots
    # shell_a = L_MAP[La_tot]
    # shell_b = L_MAP[Lb_tot]
    # return f"Cartesian ({shell_a}{shell_b}) 1-electron spin-orbit-interaction integral."

    # I think this is still faulty!
    # so1el_ints_Ls = gen_integral_exprs(
    # lambda La_tot, Lb_tot: SpinOrbitShell(
    # La_tot, Lb_tot, a, b, center_A, center_B, center_C
    # ),
    # (l_max, l_max),
    # "so1el",
    # (A_map, B_map, C_map),
    # )

    # so1el_rendered = render_py_funcs(
    # so1el_ints_Ls,
    # (a, A, b, B, C),
    # "so1el",
    # so1el_doc_func,
    # add_imports=boys_import,
    # )
    # write_file(out_dir, "so1el.py", so1el_rendered)
    # print()

    #########################
    # 2el Coulomb Integrals #
    #########################

    # def eri_doc_func(L_tots):
    # La_tot, Lb_tot, Lc_tot, Ld_tot = L_tots
    # shell_a = L_MAP[La_tot]
    # shell_b = L_MAP[Lb_tot]
    # shell_c = L_MAP[Lc_tot]
    # shell_d = L_MAP[Ld_tot]
    # return (
    # f"Cartesian [{shell_a}{shell_b}|{shell_c}{shell_d}] "
    # "2-electron electron repulsion integral."
    # )

    # I think this is still faulty!
    # eri_ints_Ls = gen_integral_exprs(
    # lambda La_tot, Lb_tot, Lc_tot, Ld_tot: ERIShell(
    # La_tot,
    # Lb_tot,
    # Lc_tot,
    # Ld_tot,
    # a,
    # b,
    # c,
    # d,
    # center_A,
    # center_B,
    # center_C,
    # center_D,
    # ),
    # # (l_max, l_max, l_max, l_max),
    # (2, 2, 2, 2),  # Stop at [dd|dd] for now
    # "eri",
    # (A_map, B_map, C_map, D_map),
    # )
    # eri_rendered = render_py_funcs(
    # eri_ints_Ls,
    # (a, A, b, B, c, C, d, D),
    # "eri",
    # eri_doc_func,
    # add_imports=boys_import,
    # )
    # write_file(out_dir, "eri.py", eri_rendered)
    # print()

    funcs = {
        "cgto": cart_gto,  # Cartesian Gaussian-type-orbital for density evaluation
        "ovlp": overlap,  # Overlap integrals
        "dpm": dipole,  # Linear moment (dipole) integrals
        "dqpm": diag_quadrupole,  # Diagonal part of the quadrupole tensor
        "qpm": quadrupole,  # Quadratic moment (quadrupole) integrals
        "kin": kinetic,  # Kinetic energy integrals
        "coul": coulomb,  # 1-electron Coulomb integrals
        # Integrals for density fitting
        "2c2e": _2center2electron,  # 2-center-2-electron integrals for DF
        "3c2e": _3center2electron,  # 3-center-2-electron integrals for density fitting
        "3c2e_sph": _3center2electron_sph,  # 3c2el for spherical transformation for DF
    }
    if len(keys) == 0:
        keys = funcs.keys()
    for key in keys:
        funcs[key]()


if __name__ == "__main__":
    run()
