#!/usr/bin/env python

"""
https://stackoverflow.com/questions/39536540/substitute-function-call-with-sympy
https://stackoverflow.com/questions/36197283/recursive-substitution-in-sympy

Functions module
http://docs.sympy.org/latest/modules/functions/index.html

CSE and code generation
http://www.sympy.org/scipy-2017-codegen-tutorial/notebooks/07-the-hard-way.html
"""

import argparse
import itertools as it
import os
from pathlib import Path
import random
import string
import sys
import textwrap
import time

from jinja2 import Template
from sympy import (
    Array,
    cse,
    exp,
    Function,
    IndexedBase,
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


def make_py_func(exprs, args=None, name=None, doc_str=""):
    if args is None:
        args = list()
    # Generate random name, if no name was supplied
    if name is None:
        name = "func_" + "".join(
            [random.choice(string.ascii_letters) for i in range(8)]
        )
    arg_strs = [arg.strip() for arg in args.split(",")]

    repls, reduced = cse(list(exprs))

    print_func = NumPyPrinter().doprint
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
            return sqrt(pi / p) * exp(-mu * X ** 2)
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
        return exprs


class Overlap3dShell(Function):
    """Whole shell of 3d overlap integrals for a given L pair
    over primitive gaussians."""

    @classmethod
    def eval(cls, La_tot, Lb_tot, a, b, A, B):
        exprs = Multipole3dShell(La_tot, Lb_tot, a, b, A, B)
        return exprs


def get_center(i):
    symbs = [Symbol(str(i) + ind) for ind in ("x", "y", "z")]
    return Array([*symbs])


def get_map(i, center_i):
    array = IndexedBase(i, shape=3)
    array_map = dict(zip(center_i, array))
    return array, array_map


def gen_integrals_exprs(int_func, L_maxs, kind):
    ranges = [range(L + 1) for L in L_maxs]
    for L_tots in it.product(*ranges):
        start = time.time()
        print(f"Generating {L_tots} {kind} ... ", end="")
        exprs = int_func(*L_tots)
        dur = time.time() - start
        print(f"finished in {dur:.2f} s!")
        yield exprs, L_tots


def render_py_funcs(exprs_Ls, args, base_name, doc_func):
    args = ", ".join((map(str, args)))

    py_funcs = list()
    for expr, L_tots in exprs_Ls:
        doc_str = doc_func(L_tots)
        doc_str += "\n\nGenerated code; DO NOT modify by hand!"
        name = base_name + "_" + "".join(str(l) for l in L_tots)
        py_func = make_py_func(expr, args=args, name=name, doc_str=doc_str)
        py_funcs.append(py_func)

    py_funcs_joined = "\n\n".join(py_funcs)
    np_tpl = Template(
        "import numpy\n\n{{ py_funcs }}",
        trim_blocks=True,
        lstrip_blocks=True,
    )
    np_rendered = np_tpl.render(py_funcs=py_funcs_joined)
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
    print(f"Wrote '{name}'.")


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("--lmax", default=L_MAX, type=int)

    return parser.parse_args()


def run():
    args = parse_args(sys.argv[1:])

    l_max = args.lmax
    out_dir = Path("gen_ints")
    try:
        os.mkdir(out_dir)
    except FileExistsError:
        pass

    # Cartesian basis function centers A and B.
    center_A = get_center("A")
    center_B = get_center("B")
    center_C = get_center("C")
    # Cartesian components (x, y, z) of the centers A and B.
    Ax, Ay, Az = center_A
    Bx, By, Bz = center_B
    # Orbital exponents a and b.
    a, b = symbols("a b")

    # Multipole moment integral center
    Cx, Cy, Cz = center_C

    # These maps will be used to convert {Ax, Ay, ...} to array quantities
    # in the generated code. This way an iterable/np.ndarray can be used as
    # function argument instead of (Ax, Ay, Az, Bx, By, Bz).
    A, A_map = get_map("A", center_A)
    B, B_map = get_map("B", center_B)
    C, C_map = get_map("C", center_C)

    def ovlp_func(La_tot, Lb_tot):
        return (
            Array(Overlap3dShell(La_tot, Lb_tot, a, b, (Ax, Ay, Az), (Bx, By, Bz)))
            .xreplace(A_map)
            .xreplace(B_map)
        )

    def ovlp_doc_func(L_tots):
        La_tot, Lb_tot = L_tots
        shell_a = L_MAP[La_tot]
        shell_b = L_MAP[Lb_tot]
        return f"Cartesian 3D ({shell_a}{shell_b}) overlap integral."

    def dipole_func(La_tot, Lb_tot, Le_tot=1):
        return (
            Array(
                Multipole3dShell(
                    La_tot,
                    Lb_tot,
                    a,
                    b,
                    (Ax, Ay, Az),
                    (Bx, By, Bz),
                    Le_tot,
                    (Cx, Cy, Cz),
                )
            )
            .xreplace(A_map)
            .xreplace(B_map)
            .xreplace(C_map)
        )

    def dipole_doc_func(L_tots):
        La_tot, Lb_tot = L_tots
        shell_a = L_MAP[La_tot]
        shell_b = L_MAP[Lb_tot]
        return f"Cartesian 3D ({shell_a}{shell_b}) dipole moment integral."

    ovlp_ints_Ls = gen_integrals_exprs(ovlp_func, (l_max, l_max), "overlap")
    ovlp_rendered = render_py_funcs(ovlp_ints_Ls, (a, A, b, B), "ovlp3d", ovlp_doc_func)
    write_py(out_dir, "ovlps3d.py", ovlp_rendered)
    print()

    dipole_ints_Ls = gen_integrals_exprs(dipole_func, (l_max, l_max), "dipole moment")
    dipole_rendered = render_py_funcs(
        dipole_ints_Ls, (a, A, b, B, C), "dipole3d", dipole_doc_func
    )
    write_py(out_dir, "dipoles3d.py", dipole_rendered)
    print()


if __name__ == "__main__":
    run()
