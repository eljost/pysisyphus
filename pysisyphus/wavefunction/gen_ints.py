#!/usr/bin/env python

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

import argparse
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
    cse,
    exp,
    Function,
    IndexedBase,
    Matrix,
    pi,
    sqrt,
    Symbol,
    symbols,
)
from sympy.codegen.ast import Assignment
from sympy.printing.numpy import NumPyPrinter

try:
    from pysisyphus.config import L_MAX
except ModuleNotFoundError:
    L_MAX = 4


L_MAP = {
    0: "s",
    1: "p",
    2: "d",
    3: "f",
    4: "g",
    5: "h",
}


def make_py_func(repls, reduced, args=None, name=None, doc_str=""):
    if args is None:
        args = list()
    # Generate random name, if no name was supplied
    if name is None:
        name = "func_" + "".join(
            [random.choice(string.ascii_letters) for i in range(8)]
        )
    arg_strs = [arg.strip() for arg in args.split(",")]

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
        S = numpy.array({{ return_val }})
        return S
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
            arg_strs=arg_strs,
        )
    ).strip()
    return rendered


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


class Multipole1d(Function):
    """1d multipole-moment integral of order 'e', between primitive 1d Gaussians
    Ga = G_i(a, r, A) and Gb = G_j(b, r, B) with Cartesian quantum number i and j,
    exponents a and b, centered at A (B). The origin of the multipole expansion is
    at C.
    """

    @classmethod
    @functools.cache
    def eval(cls, i, j, a, b, A, B, e, C):
        if i.is_negative or j.is_negative or e.is_negative:
            return 0

        p = a + b
        P = (a * A + b * B) / p

        def recur(i, j, e):
            """Simple wrapper to pass all required arguments."""
            return Multipole1d(i, j, a, b, A, B, e, C)

        def recur_rel(i, j, e, X):
            return X * recur(i, j, e) + 1 / (2 * p) * (
                i * recur(i - 1, j, e) + j * recur(i, j - 1, e) + e * recur(i, j, e - 1)
            )

        # Base case
        if i.is_zero and j.is_zero and e.is_zero:
            X = A - B
            mu = a * b / p
            return sqrt(pi / p) * exp(-mu * X**2)
        # Decrement i
        elif i.is_positive:
            X = P - A
            return recur_rel(i - 1, j, e, X)
        # Decrement j
        elif j.is_positive:
            X = P - B
            return recur_rel(i, j - 1, e, X)
        elif e.is_positive:
            X = P - C
            return recur_rel(i, j, e - 1, X)


class Multipole3d(Function):
    @classmethod
    def eval(cls, La, Lb, a, b, A, B, Le, C):
        x, y, z = [
            Multipole1d(La[i], Lb[i], a, b, A[i], B[i], Le[i], C[i]) for i in range(3)
        ]
        return x * y * z


class Multipole3dShell(Function):
    @classmethod
    def eval(cls, La_tot, Lb_tot, a, b, A, B, Le_tot=0, C=(0.0, 0.0, 0.0)):
        exprs = [
            Multipole3d(La, Lb, a, b, A, B, Le, C)
            for La, Lb, Le in shell_iter((La_tot, Lb_tot, Le_tot))
        ]
        # print(Multipole1d.eval.cache_info())
        return exprs


class Overlap1d(Function):
    @classmethod
    def eval(cls, i, j, a, b, A, B):
        return Multipole1d(i, j, a, b, A, B, 0, (0.0, 0.0, 0.0))


class Overlap3dShell(Function):
    """Whole shell of 3d overlap integrals for a given L pair
    over primitive gaussians."""

    @classmethod
    def eval(cls, La_tot, Lb_tot, a, b, A, B):
        exprs = Multipole3dShell(La_tot, Lb_tot, a, b, A, B)
        # print(Multipole1d.eval.cache_info())
        return exprs


class Kinetic1d(Function):
    @classmethod
    @functools.cache
    def eval(cls, i, j, a, b, A, B):
        if i < 0 or j < 0:
            return 0

        p = a + b
        P = (a * A + b * B) / p

        def recur(i, j):
            """Simple wrapper to pass all required arguments."""
            return Kinetic1d(i, j, a, b, A, B)

        def recur_rel(i, j, X):
            return X * recur(i, j) + 1 / (2 * p) * (
                i * recur(i - 1, j) + j * recur(i, j - 1)
            )

        def recur_ovlp(i, j):
            return Overlap1d(i, j, a, b, A, B)

        # Base case
        if i == 0 and j == 0:
            X = P - A
            return (a - 2 * a**2 * (X**2 + 1 / (2 * p))) * Overlap1d(
                i, j, a, b, A, B
            )
        # Decrement i
        elif i > 0:
            X = P - A
            # Eq. (9.3.41)
            return recur_rel(i - 1, j, X) + b / p * (
                2 * a * recur_ovlp(i, j) - i * recur_ovlp(i - 2, j)
            )
        # Decrement j
        elif j > 0:
            X = P - B
            # Eq. (9.3.41)
            return recur_rel(i, j - 1, X) + a / p * (
                2 * b * recur_ovlp(i, j) - j * recur_ovlp(i, j - 2)
            )


class Kinetic3d(Function):
    @classmethod
    def eval(cls, La, Lb, a, b, A, B):
        Tx, Ty, Tz = [Kinetic1d(La[i], Lb[i], a, b, A[i], B[i]) for i in range(3)]
        Sx, Sy, Sz = [Overlap1d(La[i], Lb[i], a, b, A[i], B[i]) for i in range(3)]
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


class Coulomb(Function):
    """Nucleus at C."""

    @classmethod
    @functools.cache
    def eval(cls, i, k, m, j, l, n, N, a, b, A, B, C):
        ang_moms = (i, k, m, j, l, n)
        if any([am < 0 for am in ang_moms]):
            return 0

        p = a + b
        P = (a * A + b * B) / p
        mu = (a * b) / p
        X_AB = A - B
        X_PC = P - C

        def recur(N, *inds):
            """Simple wrapper to pass all required arguments."""
            return Coulomb(*inds, N, a, b, A, B, C)

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

            return (
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
            )

        # Base case
        if all([am == 0 for am in ang_moms]):
            r2_PC = X_PC.dot(X_PC)
            r2_AB = X_AB.dot(X_AB)
            K = exp(-mu * r2_AB)
            return 2 * pi / p * K * boys(N, p * r2_PC)
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
            Coulomb(*La, *Lb, 0, a, b, A, B, C)
            for La, Lb in shell_iter((La_tot, Lb_tot))
        ]
        # print(Coulomb.eval.cache_info())
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


def get_map(i, center_i):
    array = IndexedBase(i, shape=3)
    array_map = dict(zip(center_i, array))
    return array, array_map


def gen_integral_exprs(
    int_func,
    L_maxs,
    kind,
    maps=None,
):
    if maps is None:
        maps = list()

    ranges = [range(L + 1) for L in L_maxs]

    for L_tots in it.product(*ranges):
        time_str = time.strftime("%H:%M:%S")
        start = time.time()
        print(f"{time_str} - Generating {L_tots} {kind} ... ", end="")
        sys.stdout.flush()
        # Generate actual expressions
        exprs = int_func(*L_tots)

        # Common subexpression elimination
        repls, reduced = cse(list(exprs))

        for i, red in enumerate(reduced):
            reduced[i] = functools.reduce(
                lambda red, map_: red.xreplace(map_), maps, red
            )

        for i, (lhs, rhs) in enumerate(repls):
            repls[i] = (
                lhs,
                functools.reduce(lambda rhs, map_: rhs.xreplace(map_), maps, rhs),
            )

        dur = time.time() - start
        print(f"finished in {dur: >8.2f} s!")
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
        py_func = make_py_func(repls, reduced, args=args, name=name, doc_str=doc_str)
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


def write_py(out_dir, name, py_rendered):
    py_name = out_dir / name
    with open(py_name, "w") as handle:
        handle.write(py_rendered)
    print(f"Wrote '{py_name}'.")


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--lmax",
        default=L_MAX,
        type=int,
        help="Generate 1e-integrals up to this maximum angular momentum.",
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

    return parser.parse_args()


def run():
    args = parse_args(sys.argv[1:])

    l_max = args.lmax
    out_dir = Path(args.out_dir if not args.write else ".")
    try:
        os.mkdir(out_dir)
    except FileExistsError:
        pass

    # Cartesian basis function centers A and B.
    center_A = get_center("A")
    center_B = get_center("B")
    center_C = get_center("C")
    center_D = get_center("D")

    # Cartesian components (x, y, z) of the centers A, B, C and D.
    # Ax, Ay, Az = center_A
    # Bx, By, Bz = center_B
    # Cx, Cy, Cz = center_C
    # Dx, Dy, Dz = center_D

    # Orbital exponents a, b, c, d.
    a, b, c, d = symbols("a b c d", real=True)

    # These maps will be used to convert {Ax, Ay, ...} to array quantities
    # in the generated code. This way an iterable/np.ndarray can be used as
    # function argument instead of (Ax, Ay, Az, Bx, By, Bz).
    A, A_map = get_map("A", center_A)
    B, B_map = get_map("B", center_B)
    C, C_map = get_map("C", center_C)
    D, D_map = get_map("D", center_D)

    #####################
    # Overlap integrals #
    #####################

    def ovlp_doc_func(L_tots):
        La_tot, Lb_tot = L_tots
        shell_a = L_MAP[La_tot]
        shell_b = L_MAP[Lb_tot]
        return f"Cartesian 3D ({shell_a}{shell_b}) overlap integral."

    ovlp_ints_Ls = gen_integral_exprs(
        lambda La_tot, Lb_tot: Overlap3dShell(La_tot, Lb_tot, a, b, center_A, center_B),
        (l_max, l_max),
        "overlap",
        (A_map, B_map),
    )
    ovlp_rendered = render_py_funcs(ovlp_ints_Ls, (a, A, b, B), "ovlp3d", ovlp_doc_func)
    write_py(out_dir, "ovlp3d.py", ovlp_rendered)
    print()

    ###########################
    # Dipole moment integrals #
    ###########################

    def dipole_doc_func(L_tots):
        La_tot, Lb_tot = L_tots
        shell_a = L_MAP[La_tot]
        shell_b = L_MAP[Lb_tot]
        return (
            f"Cartesian 3D ({shell_a}{shell_b}) dipole moment integral.\n"
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

    dipole_rendered = render_py_funcs(
        dipole_ints_Ls,
        (a, A, b, B, C),
        "dipole3d",
        dipole_doc_func,
        comment=dipole_comment,
    )
    write_py(out_dir, "dipole3d.py", dipole_rendered)
    print()

    ############################
    # Kinetic energy integrals #
    ############################

    def kinetic_doc_func(L_tots):
        La_tot, Lb_tot = L_tots
        shell_a = L_MAP[La_tot]
        shell_b = L_MAP[Lb_tot]
        return f"Cartesian 3D ({shell_a}{shell_b}) kinetic energy integral."

    kinetic_ints_Ls = gen_integral_exprs(
        lambda La_tot, Lb_tot: Kinetic3dShell(La_tot, Lb_tot, a, b, center_A, center_B),
        (l_max, l_max),
        "kinetic",
        (A_map, B_map),
    )
    kinetic_rendered = render_py_funcs(
        kinetic_ints_Ls, (a, A, b, B), "kinetic3d", kinetic_doc_func
    )
    write_py(out_dir, "kinetic3d.py", kinetic_rendered)
    print()

    #########################
    # 1el Coulomb Integrals #
    #########################

    def coulomb_doc_func(L_tots):
        La_tot, Lb_tot = L_tots
        shell_a = L_MAP[La_tot]
        shell_b = L_MAP[Lb_tot]
        return f"Cartesian ({shell_a}{shell_b}) 1-electron Coulomb integral."

    boys_import = ("from pysisyphus.wavefunction.boys import boys",)

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
    write_py(out_dir, "coulomb3d.py", coulomb_rendered)
    print()

    ####################################
    # Spin-orbit interaction integrals #
    ####################################

    def so1el_doc_func(L_tots):
        La_tot, Lb_tot = L_tots
        shell_a = L_MAP[La_tot]
        shell_b = L_MAP[Lb_tot]
        return f"Cartesian ({shell_a}{shell_b}) 1-electron spin-orbit-interaction integral."

    so1el_ints_Ls = gen_integral_exprs(
        lambda La_tot, Lb_tot: SpinOrbitShell(
            La_tot, Lb_tot, a, b, center_A, center_B, center_C
        ),
        (l_max, l_max),
        "so1el",
        (A_map, B_map, C_map),
    )

    so1el_rendered = render_py_funcs(
        so1el_ints_Ls,
        (a, A, b, B, C),
        "so1el",
        so1el_doc_func,
        add_imports=boys_import,
    )
    write_py(out_dir, "so1el.py", so1el_rendered)
    print()

    #########################
    # 2el Coulomb Integrals #
    #########################

    def eri_doc_func(L_tots):
        La_tot, Lb_tot, Lc_tot, Ld_tot = L_tots
        shell_a = L_MAP[La_tot]
        shell_b = L_MAP[Lb_tot]
        shell_c = L_MAP[Lc_tot]
        shell_d = L_MAP[Ld_tot]
        return (
            f"Cartesian [{shell_a}{shell_b}|{shell_c}{shell_d}] "
            "2-electron electron repulsion integral."
        )

    eri_ints_Ls = gen_integral_exprs(
        lambda La_tot, Lb_tot, Lc_tot, Ld_tot: ERIShell(
            La_tot,
            Lb_tot,
            Lc_tot,
            Ld_tot,
            a,
            b,
            c,
            d,
            center_A,
            center_B,
            center_C,
            center_D,
        ),
        # (l_max, l_max, l_max, l_max),
        (2, 2, 2, 2),  # Stop at [dd|dd] for now
        "eri",
        (A_map, B_map, C_map, D_map),
    )
    eri_rendered = render_py_funcs(
        eri_ints_Ls,
        (a, A, b, B, c, C, d, D),
        "eri",
        eri_doc_func,
        add_imports=boys_import,
    )
    write_py(out_dir, "eri.py", eri_rendered)
    print()


if __name__ == "__main__":
    run()
