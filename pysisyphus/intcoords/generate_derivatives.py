#!/usr/bin/env python3
# Johannes Steinmetzer 2019
# As describend in:
#   V. Bakken, T. Helgaker, J. Chem. Phys., 117, 20, 2002
# [1] https://aip.scitation.org/doi/abs/10.1063/1.1515483
# [2] https://doi.org/10.1002/(SICI)1096-987X(19960115)17:1<49::AID-JCC5>3.0.CO;2-0

from collections import namedtuple
import random
import string
import textwrap

from jinja2 import Template
import numpy as np
import sympy as sym
from sympy import cse
from sympy.codegen.ast import Assignment
from sympy.printing.pycode import pycode
from sympy.vector import CoordSys3D


Derivs = namedtuple("Derivs",
                    "d1 d2 f1 f2",
)


def make_py_func(exprs, args=None, name=None, comment=""):
    if args is None:
        args = list()
    if name is None:
        name = "func_" + "".join([random.choice(string.ascii_letters)
                                  for i in range(8)])

    repls, reduced = cse(list(exprs))

    assignments = [Assignment(lhs, rhs) for lhs, rhs in repls]
    py_lines = [pycode(as_) for as_ in assignments]
    return_val = pycode(reduced)

    tpl = Template("""
    def {{ name }}({{ args }}):
        {% if comment %}
        \"\"\"{{ comment }}\"\"\"
        {% endif %}

        {% for line in py_lines %}
        {{ line }}
        {% endfor %}

        return np.array({{ return_val }})
    """, trim_blocks=True, lstrip_blocks=True)

    rendered = textwrap.dedent(
                    tpl.render(
                        name=name,
                        args=args,
                        py_lines=py_lines,
                        return_val=return_val,
                        comment=comment,
                    )
    ).strip()
    return rendered


def make_deriv_funcs(base_expr, dx, args, names, comments):
    d1_name, d2_name = names
    d1_comment, d2_comment = comments
    # First derivative
    deriv1 = sym.derive_by_array(base_expr, dx)
    deriv1_func = make_py_func(deriv1, args=args, name=d1_name,
                               comment=d1_comment)
    # Second derivative
    deriv2 = sym.derive_by_array(deriv1, dx)
    deriv2_func = make_py_func(deriv2, args=args, name=d2_name,
                              comment=d2_comment)
    return Derivs(
        d1=deriv1,
        d2=deriv2,
        f1=deriv1_func,
        f2=deriv2_func,
    )


def generate_wilson():
    m0, m1, m2, n0, n1, n2, o0, o1, o2, p0, p1, p2 = sym.symbols("m:3 n:3 o:3 p:3")

    # Coordinate system
    Sys = CoordSys3D("Sys")
    M = Sys.origin.locate_new("M", m0*Sys.i + m1*Sys.j + m2*Sys.k)
    N = Sys.origin.locate_new("N", n0*Sys.i + n1*Sys.j + n2*Sys.k)
    O = Sys.origin.locate_new("O", o0*Sys.i + o1*Sys.j + o2*Sys.k)
    P = Sys.origin.locate_new("P", p0*Sys.i + p1*Sys.j + p2*Sys.k)

    # Bond/Stretch
    U = M.position_wrt(N)
    q_b = U.magnitude()
    dx_b = (m0, m1, m2, n0, n1, n2)
    args_b = "m0, m1, m2, n0, n1, n2"
    derivs_b = make_deriv_funcs(q_b, dx_b, args_b,
                                ("dq_b", "d2q_b"),
                                ("Stretch, first derivative wrt. cartesians",
                                 "Stretch, 2nd derivative wrt. cartesians",)
    )
    print(derivs_b.f1)
    print(derivs_b.f2)

    # Bend/Angle
    U = M.position_wrt(O)
    V = N.position_wrt(O)
    q_a = sym.acos(U.dot(V) / (U.magnitude() * V.magnitude()))
    dx_a = (m0, m1, m2, o0, o1, o2, n0, n1, n2)
    args_a = "m0, m1, m2, o0, o1, o2, n0, n1, n2"
    derivs_a = make_deriv_funcs(q_a, dx_a, args_a,
                                ("dq_a", "d2q_a"),
                                ("Bend, first derivative wrt. cartesians",
                                 "Bend, 2nd derivative wrt. cartesians",)
    )
    print(derivs_a.f1)
    print(derivs_a.f2)

    # Dihedral/Torsion
    U = M.position_wrt(O)
    V = N.position_wrt(P)
    W = P.position_wrt(O)
    U_ = U.normalize()
    W_ = W.normalize()
    V_ = V.normalize()
    phi_u = sym.acos(U_.dot(W_))
    phi_v = sym.acos(-W_.dot(V_))
    q_d = sym.acos(U_.cross(W_).dot(V_.cross(W_))/(sym.sin(phi_u)*sym.sin(phi_v)))
    dx_d = (m0, m1, m2, o0, o1, o2, p0, p1, p2, n0, n1, n2)
    args_d = "m0, m1, m2, o0, o1, o2, p0, p1, p2, n0, n1, n2"
    derivs_d = make_deriv_funcs(q_d, dx_d, args_d,
                                ("dq_d", "d2q_d"),
                                ("Torsion, first derivative wrt. cartesians",
                                 "Torsion, 2nd derivative wrt. cartesians",)
    )
    print(derivs_d.f1)
    print(derivs_d.f2)

    out_fn = "derivatives.py"
    with open(out_fn, "w") as handle:
        handle.write("from math import sqrt\n\nimport numpy as np\n\n\n")
        for d in (derivs_b, derivs_a, derivs_d):
            handle.write(d.f1 + "\n\n\n")
            handle.write(d.f2 + "\n\n\n")


if __name__ == "__main__":
    generate_wilson()
