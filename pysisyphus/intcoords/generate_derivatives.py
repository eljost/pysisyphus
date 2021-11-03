#!/usr/bin/env python3

# Johannes Steinmetzer 2020
# As describend in:
#   V. Bakken, T. Helgaker, J. Chem. Phys., 117, 20, 2002
# [1] https://aip.scitation.org/doi/abs/10.1063/1.1515483
# [2] https://doi.org/10.1002/(SICI)1096-987X(19960115)17:1<49::AID-JCC5>3.0.CO;2-0

from collections import namedtuple
import itertools as it
import random
import string
import textwrap

from jinja2 import Template
import sympy as sym
from sympy import cse
from sympy.codegen.ast import Assignment
from sympy.printing.pycode import pycode, MpmathPrinter
from sympy.vector import CoordSys3D
from sympy.tensor.array.dense_ndim_array import ImmutableDenseNDimArray


FuncResult = namedtuple("FuncResult", "d0 d1 d2 f0 f1 f2")


def make_py_func(exprs, args=None, name=None, comment="", use_mpmath=False):
    if args is None:
        args = list()
    if name is None:
        name = "func_" + "".join(
            [random.choice(string.ascii_letters) for i in range(8)]
        )
    arg_strs = [arg.strip() for arg in args.split(",")]

    is_scalar = not isinstance(exprs, ImmutableDenseNDimArray)
    if is_scalar:
        repls, reduced = cse(exprs)
        reduced = reduced[0]
    else:
        if len(exprs.shape) == 2:
            exprs = it.chain(*exprs)
        repls, reduced = cse(list(exprs))

    print_func = pycode
    if use_mpmath:
        print_func = MpmathPrinter().doprint

    assignments = [Assignment(lhs, rhs) for lhs, rhs in repls]
    py_lines = [print_func(as_) for as_ in assignments]
    return_val = print_func(reduced)

    tpl = Template(
        """
    def {{ name }}({{ args }}):
        {% if comment %}
        \"\"\"{{ comment }}\"\"\"
        {% endif %}

        {% if use_mpmath %}
        {% for arg in arg_strs %}
        {{ arg }} = mpmath.mpf({{ arg }})
        {% endfor %}
        {% endif %}

        {% for line in py_lines %}
        {{ line }}
        {% endfor %}

        {% if is_scalar %}
        return {{ return_val }}
        {% else %}
        return np.array({{ return_val }}, dtype=np.float64)
        {% endif %}
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
            comment=comment,
            is_scalar=is_scalar,
            use_mpmath=use_mpmath,
            arg_strs=arg_strs,
        )
    ).strip()
    return rendered


# def make_deriv_funcs(base_expr, dx, args, names, comments):
def make_deriv_funcs(base_expr, dx, args, names, comment, use_mpmath=True):
    q_name, d1_name, d2_name = names
    q_comment, d1_comment, d2_comment = [
        comment + add
        for add in [
            "",
            ", first derivative wrt. Cartesians",
            ", 2nd derivative wrt. Cartesians",
        ]
    ]

    # Actual function
    print(f"\tmpmath: {use_mpmath}")
    print("\tFunction")
    q_func = make_py_func(
        base_expr, args=args, name=q_name, comment=q_comment, use_mpmath=use_mpmath
    )

    # First derivative
    print("\t1st derivative")
    deriv1 = sym.derive_by_array(base_expr, dx)
    deriv1_func = make_py_func(
        deriv1, args=args, name=d1_name, comment=d1_comment, use_mpmath=use_mpmath
    )

    # Second derivative
    print("\t2nd derivative")
    deriv2 = sym.derive_by_array(deriv1, dx)
    deriv2_func = make_py_func(
        deriv2, args=args, name=d2_name, comment=d2_comment, use_mpmath=use_mpmath
    )

    return FuncResult(
        # Expressions
        d0=base_expr,
        d1=deriv1,
        d2=deriv2,
        # Functions
        f0=q_func,
        f1=deriv1_func,
        f2=deriv2_func,
    )


def generate_wilson(generate=None, out_fn="derivatives.py", use_mpmath=False):
    m0, m1, m2, n0, n1, n2, o0, o1, o2, p0, p1, p2 = sym.symbols("m:3 n:3 o:3 p:3")

    # Coordinate system
    Sys = CoordSys3D("Sys")
    M = Sys.origin.locate_new("M", m0 * Sys.i + m1 * Sys.j + m2 * Sys.k)
    N = Sys.origin.locate_new("N", n0 * Sys.i + n1 * Sys.j + n2 * Sys.k)
    O = Sys.origin.locate_new("O", o0 * Sys.i + o1 * Sys.j + o2 * Sys.k)
    P = Sys.origin.locate_new("P", p0 * Sys.i + p1 * Sys.j + p2 * Sys.k)

    def bond():
        # Bond/Stretch
        U = M.position_wrt(N)
        q_b = U.magnitude()
        dx_b = (m0, m1, m2, n0, n1, n2)
        args_b = "m0, m1, m2, n0, n1, n2"
        func_result_b = make_deriv_funcs(
            q_b,
            dx_b,
            args_b,
            ("q_b", "dq_b", "d2q_b"),
            "Stretch",
            use_mpmath=use_mpmath,
        )
        return func_result_b

    def bend():
        # Bend/Angle
        U = M.position_wrt(O)
        V = N.position_wrt(O)
        q_a = sym.acos(U.dot(V) / (U.magnitude() * V.magnitude()))
        dx_a = (m0, m1, m2, o0, o1, o2, n0, n1, n2)
        args_a = "m0, m1, m2, o0, o1, o2, n0, n1, n2"
        func_result_a = make_deriv_funcs(
            q_a,
            dx_a,
            args_a,
            ("q_a", "dq_a", "d2q_a"),
            "Bend",
            use_mpmath=use_mpmath,
        )
        return func_result_a

    def bend2():
        # atan2 based Bend/Angle
        U = M.position_wrt(O)
        V = N.position_wrt(O)
        q_a2 = sym.atan2(U.cross(V).magnitude(), U.dot(V))
        dx_a2 = (m0, m1, m2, o0, o1, o2, n0, n1, n2)
        args_a2 = "m0, m1, m2, o0, o1, o2, n0, n1, n2"
        func_result_a = make_deriv_funcs(
            q_a2,
            dx_a2,
            args_a2,
            ("q_a2", "dq_a2", "d2q_a2"),
            "Bend2",
            use_mpmath=use_mpmath,
        )
        return func_result_a

    def dihedral():
        # Dihedral/Torsion
        U = M.position_wrt(O)
        W = P.position_wrt(O)
        V = N.position_wrt(P)
        U_ = U.normalize()
        W_ = W.normalize()
        V_ = V.normalize()
        phi_u = sym.acos(U_.dot(W_))
        phi_v = sym.acos(-W_.dot(V_))
        q_d = sym.acos(
            U_.cross(W_).dot(V_.cross(W_)) / (sym.sin(phi_u) * sym.sin(phi_v))
        )
        dx_d = (m0, m1, m2, o0, o1, o2, p0, p1, p2, n0, n1, n2)
        args_d = "m0, m1, m2, o0, o1, o2, p0, p1, p2, n0, n1, n2"
        func_result_d = make_deriv_funcs(
            q_d,
            dx_d,
            args_d,
            ("q_d", "dq_d", "d2q_d"),
            "Torsion",
            use_mpmath=use_mpmath,
        )
        return func_result_d

    def dihedral2():
        # atan2 based Dihedral/Torsion
        U1 = O.position_wrt(M)
        U2 = P.position_wrt(O)
        U3 = N.position_wrt(P)
        cross_U2U3 = U2.cross(U3)
        q_d2 = sym.atan2((U2.magnitude() * U1).dot(cross_U2U3), cross_U2U3.dot(U1.cross(U2)))
        dx_d2 = (m0, m1, m2, o0, o1, o2, p0, p1, p2, n0, n1, n2)
        args_d2 = "m0, m1, m2, o0, o1, o2, p0, p1, p2, n0, n1, n2"
        func_result_d = make_deriv_funcs(
            q_d2,
            dx_d2,
            args_d2,
            ("q_d2", "dq_d2", "d2q_d2"),
            "Torsion2",
            use_mpmath=use_mpmath,
        )
        return func_result_d

    def linear_bend():
        # Linear Bend
        U = M.position_wrt(O)
        V = N.position_wrt(O)
        W = P.position_wrt(Sys)
        q_lb = W.dot(U.cross(V)) / (U.magnitude() * V.magnitude())
        dx_lb = (m0, m1, m2, o0, o1, o2, n0, n1, n2)
        # Additional args, as we also supply an orthogonal direction
        args_lb = "m0, m1, m2, o0, o1, o2, n0, n1, n2, p0, p1, p2"
        func_result_lb = make_deriv_funcs(
            q_lb,
            dx_lb,
            args_lb,
            ("q_lb", "dq_lb", "d2q_lb"),
            "Linear Bend",
            use_mpmath=use_mpmath,
        )
        return func_result_lb

    def out_of_plane():
        U = M.position_wrt(P)
        V = N.position_wrt(P)
        W = O.position_wrt(P)

        U_ = U.normalize()
        W_ = W.normalize()
        V_ = V.normalize()

        Z = U_.cross(V_) + V_.cross(W_) + W_.cross(U_)
        Z_ = Z.normalize()

        q_oop = Z_.dot(U_)

        dx_oop = (m0, m1, m2, n0, n1, n2, o0, o1, o2, p0, p1, p2)
        args_oop = "m0, m1, m2, n0, n1, n2, o0, o1, o2, p0, p1, p2"
        func_result_oop = make_deriv_funcs(
            q_oop,
            dx_oop,
            args_oop,
            ("q_oop", "dq_oop", "d2q_oop"),
            "OutOfPlane",
            use_mpmath=use_mpmath,
        )
        return func_result_oop

    def linear_displacement():
        U = M.position_wrt(O)
        V = N.position_wrt(O)
        W = N.position_wrt(M)
        U_ = U.normalize()
        V_ = V.normalize()
        W_ = W.normalize()

        # Vector for cross product. For the complement X should correspond
        # to the first orthogonal direction (X = W_.cross(X)).
        X = P.position_wrt(Sys)
        # Orthogonal direction
        Y = W_.cross(X)
        Y_ = Y.normalize()

        q_ld = Y_.dot(U_) + Y_.dot(V_)
        dx_ld = (m0, m1, m2, o0, o1, o2, n0, n1, n2)
        # Additional args, as we also supply an orthogonal direction
        args_ld = "m0, m1, m2, o0, o1, o2, n0, n1, n2, p0, p1, p2"
        func_result_ld = make_deriv_funcs(
            q_ld,
            dx_ld,
            args_ld,
            ("q_ld", "dq_ld", "d2q_ld"),
            "Linear Displacement",
            use_mpmath=use_mpmath,
        )
        return func_result_ld

    if generate is None:
        generate = (
            "bond",
            "bend",
            "bend2",
            "dihedral",
            "dihedral2",
            "linear_bend",
            "out_of_plane",
            "linear_displacement",
        )

    avail_funcs = {
        "bond": bond,
        "bend": bend,
        "bend2": bend2,
        "dihedral": dihedral,
        "dihedral2": dihedral2,
        "linear_bend": linear_bend,
        "out_of_plane": out_of_plane,
        "linear_displacement": linear_displacement,
    }
    funcs = [avail_funcs[key] for key in generate]
    func_results = list()
    for name, func in zip(generate, funcs):
        print(f"Generating expressions for: '{name}'")
        func_res = func()
        func_results.append(func_res)

    import_str = "import mpmath" if use_mpmath else "import math"
    if out_fn:
        with open(out_fn, "w") as handle:
            handle.write(f"{import_str}\n\nimport numpy as np\n\n\n")
            for fr in func_results:
                handle.write(fr.f0 + "\n\n\n")
                handle.write(fr.f1 + "\n\n\n")
                handle.write(fr.f2 + "\n\n\n")
    print(f"Wrote generated code to '{out_fn}'")

    return func_results


if __name__ == "__main__":
    generate_wilson(out_fn="derivatives.py", use_mpmath=False)
    # print()
    generate_wilson(out_fn="mp_derivatives.py", use_mpmath=True)
    # generate_wilson(out_fn="lindisp.py", use_mpmath=False)
