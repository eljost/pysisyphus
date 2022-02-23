from pprint import pprint as pp_

import pyparsing as pp
import numpy as np

from pysisyphus.constants import ANG2BOHR
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import file_or_str


@file_or_str(".mol2")
def parse_mol2(text):
    molecule = (
        pp.CaselessLiteral("@<TRIPOS>MOLECULE")
        + pp.Word(pp.printables).set_results_name("mol_name")
        + pp.common.integer.set_results_name("num_atoms")
        + pp.Optional(pp.common.integer.set_results_name("num_bonds"))
        + pp.Optional(pp.common.integer.set_results_name("num_subst"))
        + pp.Optional(pp.common.integer.set_results_name("num_feat"))
        + pp.Optional(pp.common.integer.set_results_name("num_sets"))
        + pp.oneOf(
            ("SMALL", "BIOPOLYMER", "PROTEIN", "NUCLEIC_ACID", "SACCHARIDE")
        ).set_results_name("mol_type")
        + pp.oneOf(
            ("NO_CHARGES", "GASTEIGER", "MMFF94_CHARGES", "USER_CHARGES", "HUCKEL")
        ).set_results_name("charge_type")
        + pp.Optional(pp.Word(pp.printables).set_results_name("status_bits"))
        # Only spaces may not be enough for mol comment
        + pp.Optional(pp.Word(pp.printables + " ")).set_results_name("mol_comment")
    )

    def get_line_word(*args):
        return pp.Word(*args).setWhitespaceChars(" \t")

    atom_data_line = pp.Group(
        pp.common.integer.set_results_name("atom_id")
        + pp.Word(pp.alphas).set_results_name("atom_name")
        + pp.Group(pp.common.real + pp.common.real + pp.common.real).set_results_name(
            "xyz"
        )
        + pp.Word(pp.printables).set_results_name("atom_type")
        + pp.Optional(pp.common.integer.set_results_name("subst_id"))
        + pp.Optional(pp.Word(pp.printables).set_results_name("subst_name"))
        + pp.Optional(pp.common.real.set_results_name("charge"))
        + pp.Optional(get_line_word(pp.printables).set_results_name("status_bit"))
    )
    atom = pp.CaselessLiteral("@<TRIPOS>ATOM") + pp.OneOrMore(
        atom_data_line
    ).set_results_name("atoms_xyzs")

    parser = molecule + atom
    parser.ignore(pp.helpers.python_style_comment)

    result = parser.parseString(text)
    as_dict = result.asDict()
    atoms = list()
    coords = np.zeros((as_dict["num_atoms"], 3), dtype=float)
    atoms_xyzs = as_dict["atoms_xyzs"]
    for i, line in enumerate(as_dict["atoms_xyzs"]):
        assert i + 1 == line["atom_id"]
        atoms.append(line["atom_name"])
        coords[i] = line["xyz"]

    coords *= ANG2BOHR

    return atoms, coords


def geom_from_mol2(text, **kwargs):
    atoms, coords = parse_mol2(text)
    geom = Geometry(atoms, coords, **kwargs)
    return geom
