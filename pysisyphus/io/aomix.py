import numpy as np
import pyparsing as pp

from pysisyphus.integrals import get_l, Shell, Shells


def parse_aomix(text):
    real = pp.common.real
    sreal = pp.common.sci_real
    int_ = pp.common.integer
    header = pp.CaselessLiteral("[AOMix Format]")
    title = pp.CaselessKeyword("[Title]")
    energy = pp.Suppress(pp.CaselessLiteral("[SCF Energy / Hartree]")) + sreal
    xyz = pp.Group(sreal + sreal + sreal)
    atom_line = pp.Group(
        pp.Word(pp.alphas).set_results_name("atom")
        + int_.set_results_name("id")
        + int_.set_results_name("atom_num")
        + xyz.set_results_name("xyz")
    )
    atoms = (
        pp.Suppress(pp.CaselessLiteral("[Atoms]"))
        + pp.Optional(pp.CaselessLiteral("AU")).set_results_name("unit")
        + pp.OneOrMore(atom_line).set_results_name("atom_lines")
    )
    prim_gto = pp.Group(sreal + sreal)
    contr_gto = pp.Group(
        pp.Word(pp.alphas).set_results_name("l")
        + int_.set_results_name("contr_depth")
        + real
        + pp.OneOrMore(prim_gto).set_results_name("prim_gtos")
    )
    center_contr_gto = pp.Group(
        int_.set_results_name("center")
        + int_
        + pp.OneOrMore(contr_gto).set_results_name("contr_gtos")
    )
    gtos = pp.CaselessLiteral("[GTO]") + pp.OneOrMore(
        center_contr_gto
    ).set_results_name("center_gtos")

    parser = (
        header
        + pp.Optional(title).set_results_name("title")
        + pp.Optional(energy).set_results_name("energy")
        + atoms
        + gtos
    )

    result = parser.parseString(text)
    as_dict = result.asDict()

    atom_lines = as_dict["atom_lines"]
    center_gtos = as_dict["center_gtos"]
    _shells = list()
    for atom_line, cgtos in zip(atom_lines, center_gtos):
        center_ind = cgtos["center"]
        assert atom_line["id"] == center_ind
        center = np.array(atom_line["xyz"])
        for cgto in cgtos["contr_gtos"]:
            exps, coeffs = zip(*cgto["prim_gtos"])
            exps = np.array(exps)
            coeffs = np.array(coeffs)
            L = get_l(cgto["l"])
            shell = Shell(
                L=L,
                center=center,
                coeffs=coeffs,
                exps=exps,
                center_ind=center_ind,
            )
            _shells.append(shell)
    shells = Shells(_shells)
    return shells
