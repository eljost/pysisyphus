import jinja2
import pyparsing as pp
import numpy as np

from pysisyphus.constants import ANG2BOHR
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import file_or_str


def optional_on_line(pattern):
    return pp.Optional(~pp.LineEnd() + pattern)


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
        + pp.Word(pp.alphanums).set_results_name("atom_name")
        + pp.Group(pp.common.real + pp.common.real + pp.common.real).set_results_name(
            "xyz"
        )
        + pp.Word(pp.printables).set_results_name("atom_type")
        # TODO: fix this parser, as it breaks when the optional values are not provided.
        # If subst_id is missing the atom_id is on the next line is mistaken for it, so
        # a fix would be to restrict parsing of this token to one line.
        + optional_on_line(pp.common.integer.set_results_name("subst_id"))
        + optional_on_line(pp.Word(pp.printables).set_results_name("subst_name"))
        + optional_on_line(pp.common.real.set_results_name("charge"))
        + optional_on_line(get_line_word(pp.printables).set_results_name("status_bit"))
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

    return result


@file_or_str(".mol2")
def atoms_coords_from_mol2(text):
    result = parse_mol2(text)
    as_dict = result.asDict()

    atoms = list()
    coords = np.zeros((as_dict["num_atoms"], 3), dtype=float)
    atoms_xyzs = as_dict["atoms_xyzs"]
    for i, line in enumerate(atoms_xyzs):
        assert i + 1 == line["atom_id"]
        atom = line["atom_type"].split(".")[0]
        atoms.append(atom)
        coords[i] = line["xyz"]
    assert len(atoms) == len(coords)

    coords *= ANG2BOHR

    return atoms, coords


def geom_from_mol2(text, **kwargs):
    atoms, coords = atoms_coords_from_mol2(text)
    geom = Geometry(atoms, coords, **kwargs)
    return geom


MOL2_DICT_TPL = jinja2.Template(
    """@<TRIPOS>MOLECULE
{{ mol_name }}
{{ num_atoms }} {{ num_bonds}} {{ num_subst }} {{ num_feat }} {{ num_sets }}
{{ mol_type }}
{{ charge_type }}

@<TRIPOS>ATOM
{%- for ax in atoms_xyzs %}
{{ ax["atom_id"] }} {{ ax["atom_name"] }} {{ render_xyz(ax["xyz"])}} {{ ax["atom_type"] }}  {{ ax["subst_id"] }} {{ ax["subst_name"] }} {{ ax["charge"]}}
{%- endfor %}
@<TRIPOS>BOND
{%- for bond in bonds %}
 {{ bond["bond_id"] }} {{ bond["origin_atom_id"]}} {{ bond["target_atom_id"] }} {{ bond["bond_type"]}}
{%- endfor %}

"""
)


def render_xyz(xyz):
    fmt = " >12.6f"
    return " ".join(map(lambda _: f"{_:{fmt}}", xyz))


def dict_to_mol2_string(as_dict):
    d = as_dict
    rendered = MOL2_DICT_TPL.render(
        mol_name=d["mol_name"],
        num_bonds=d["num_bonds"],
        num_subst=d["num_subst"],
        num_feat=d["num_feat"],
        num_sets=d["num_sets"],
        num_atoms=d["num_atoms"],
        mol_type=d["mol_type"],
        charge_type=d["charge_type"],
        atoms_xyzs=d["atoms_xyzs"],
        bonds=d["bond"],
        render_xyz=render_xyz,
    )
    return rendered
