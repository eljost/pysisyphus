
import pyparsing as pp
import numpy as np

from pysisyphus.constants import ANG2BOHR
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import file_or_str


@file_or_str(".mol2")
def parse_mol2(text):
    def get_line_word(*args):
        return pp.Word(*args).setWhitespaceChars(" \t")

    new_record = pp.CaselessLiteral("@<TRIPOS>")

    # <TRIPOS>MOLECULE
    molecule = (
        pp.CaselessLiteral("@<TRIPOS>MOLECULE")
        + pp.LineEnd()
        # comment line/molecule name
        + pp.OneOrMore(get_line_word(pp.printables)).set_results_name("mol_name")
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
        + pp.Optional(
            ~new_record + pp.Word(pp.printables).set_results_name("status_bits")
        )
        # Only spaces may not be enough for mol comment
        + pp.Optional(~new_record + pp.Word(pp.printables + " ")).set_results_name(
            "mol_comment"
        )
    )

    # <TRIPOS>ATOM
    atom_data_line = pp.Group(
        ~new_record
        + pp.common.integer.set_results_name("atom_id")
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

    # <TRIPOS>BOND
    bond_data_line = pp.Group(
        ~new_record
        + pp.common.integer.set_results_name("bond_id")
        + pp.common.integer.set_results_name("origin_atom_id")
        + pp.common.integer.set_results_name("target_atom_id")
        + pp.Word(pp.printables).set_results_name("bond_type")
    )
    bond = pp.CaselessLiteral("@<TRIPOS>BOND") + pp.OneOrMore(
        bond_data_line
    ).set_results_name("bond")

    parser = molecule + atom + pp.Optional(bond)
    parser.ignore(pp.helpers.python_style_comment)

    result = parser.parseString(text)
    as_dict = result.asDict()
    # from pprint import pprint as pp_
    # pp_(as_dict)

    atoms = list()
    coords = np.zeros((as_dict["num_atoms"], 3), dtype=float)
    atoms_xyzs = as_dict["atoms_xyzs"]
    for i, line in enumerate(atoms_xyzs):
        assert i + 1 == line["atom_id"]
        atoms.append(line["atom_name"])
        coords[i] = line["xyz"]

    coords *= ANG2BOHR

    return atoms, coords


def geom_from_mol2(text, **kwargs):
    atoms, coords = parse_mol2(text)
    geom = Geometry(atoms, coords, **kwargs)
    return geom
