import re

import numpy as np
import pyparsing as pp

from pysisyphus.Geometry import Geometry
from pysisyphus.io.xyz import parse_xyz
from pysisyphus.helpers_pure import file_or_str


@file_or_str(".molden", ".input")
def parse_molden_atoms_coords(molden_str):
    _, geometries = re.split(r"\[GEOMETRIES\] \(XYZ\)", molden_str)
    atoms_coords, comments = parse_xyz(geometries, with_comment=True)

    return atoms_coords, comments


def geoms_from_molden(fn, **kwargs):
    atoms_coords, comments = parse_molden_atoms_coords(fn)
    geoms = [
        Geometry(atoms, coords.flatten(), comment=comment, **kwargs)
        for (atoms, coords), comment in zip(atoms_coords, comments)
    ]
    return geoms


@file_or_str(".molden", ".input")
def parse_molden(text, with_mos=True):
    def get_line_word(*args):
        return pp.Word(*args).setWhitespaceChars(" \t")

    int_ = pp.common.integer
    real = pp.common.real
    sci_real = pp.common.sci_real

    molden_format = pp.CaselessLiteral("[Molden Format]")
    atoms_header = pp.CaselessLiteral("[Atoms]")
    title = pp.CaselessLiteral("[Title]") + pp.Optional(
        pp.ZeroOrMore(~atoms_header + pp.Word(pp.printables))
    )

    # [Atoms]
    # element_name number atomic_number x y z
    unit = pp.one_of("AU Angs", caseless=True)
    atom_symbol = pp.Word(pp.alphas)
    xyz = pp.Group(sci_real + sci_real + sci_real)
    atom_line = pp.Group(
        atom_symbol.set_results_name("symbol")
        + int_.set_results_name("number")
        + int_.set_results_name("atomic_number")
        + xyz.set_results_name("xyz")
    )
    atoms = (
        atoms_header
        + pp.Optional(unit).set_results_name("unit")
        + pp.OneOrMore(atom_line).set_results_name("atoms")
    )

    def pgto_exps_coeffs(s, loc, toks):
        data = toks.asList()
        exps, coeffs = np.array(data, dtype=float).reshape(-1, 2).T
        pr = pp.ParseResults.from_dict(
            {
                "exponents": exps,
                "coeffs": coeffs,
            }
        )
        # Somehow this function will result in exponents and coeffs being both
        # stored in lists of length 1.
        return pr

    # [GTO]
    ang_mom = pp.one_of("s p sp d f g h i", caseless=True)
    pgto = sci_real + sci_real
    shell = pp.Group(
        ang_mom.set_results_name("ang_mom")
        + int_.set_results_name("contr_depth")
        + real.suppress()  # Usually 1.0
        + pp.OneOrMore(pgto).set_parse_action(pgto_exps_coeffs)
    )
    atom_gtos = pp.Group(
        int_.set_results_name("number")  # Center
        + pp.Literal("0")
        + pp.OneOrMore(shell).set_results_name("shells")
    )
    gto = pp.CaselessLiteral("[GTO]") + pp.Group(
        pp.OneOrMore(atom_gtos)
    ).set_results_name("gto")

    # [5D] [7F] [9G] etc.
    ang_mom_flags = pp.ZeroOrMore(pp.one_of("[5D] [7F] [9G]", caseless=True))

    # [MO]
    sym_label = pp.Word(pp.printables)
    spin = pp.one_of("Alpha Beta", caseless=True)
    mo_coeff = pp.Suppress(int_) + real
    mo = pp.Group(
        pp.CaselessLiteral("Sym=").suppress()
        + sym_label.set_results_name("sym")
        + pp.CaselessLiteral("Ene=").suppress()
        + sci_real.set_results_name("ene")
        + pp.CaselessLiteral("Spin=").suppress()
        + spin.set_results_name("spin")
        + pp.CaselessLiteral("Occup=").suppress()
        + real.set_results_name("occup")
        + pp.Group(pp.OneOrMore(mo_coeff)).set_results_name("coeffs")
    )
    mos = pp.CaselessLiteral("[MO]") + pp.OneOrMore(mo).set_results_name("mos")

    parser = molden_format + title + atoms + gto + ang_mom_flags
    if with_mos:
        parser += mos

    res = parser.parse_string(text)
    return res.asDict()
