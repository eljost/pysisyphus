from collections import OrderedDict
import re
import struct
import textwrap

from jinja2 import Template
import numpy as np

from pysisyphus.elem_data import KNOWN_ATOMS
from pysisyphus.constants import ANG2BOHR
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import chunks
from pysisyphus.intcoords.setup_fast import find_bonds


AMINOS = (
    "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU"
    "MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()
)


def get_parser(widths):
    fmt = " ".join("{}{}".format(width, "s") for width in widths)
    fieldstruct = struct.Struct(fmt)
    parse = lambda line: tuple(s.decode() for s in fieldstruct.unpack(line.encode()))
    return parse


STRIP_RE = re.compile(r"[\d\s]*")
NAME_MAP = {
    "hh": "H",
    "he": "H",
}


def parse_atom_name(name):
    org_name = name
    assert len(name) == 4
    name = name[:2]
    stripped = STRIP_RE.sub("", name).lower()
    try:
        mapped = NAME_MAP[stripped]
    except KeyError:
        assert stripped in KNOWN_ATOMS, f"Could not parse atom name '{org_name}'"
        mapped = stripped
    return mapped.capitalize()


def parse_pdb(fn):
    with open(fn) as handle:
        atm_lines = list()
        for line in handle:
            if line.startswith("HETATM") or line.startswith("ATOM"):
                atm_lines.append(line.strip())
                continue
            elif line.startswith("MASTER"):
                master_line = line.strip()

    """
        0  Record name
        1  serial
        2  name
        3  altLoc
        4  ResName
        5  chainID
        6  resSeq
        7  iCode
        8  x
        9  y
        10 z
        11 occupancy
        12 tempFactor
        13 element
        (14 charge), not parsed

        See https://stackoverflow.com/questions/4914008
    """
    atm_widths = (6, 6, 4, 1, 4, 1, 4, 1, 11, 8, 8, 6, 6, 12)
    atm_parse = get_parser(atm_widths)

    """
        0  Record name
        1  Num. of REMARK
        2  "0"
        3  Num. of HET
        4  Num. of HELIX
        5  Num. of SHEET
        6  deprecated
        7  Num. of SITE
        8  numXform
        9  Num. of atomic coordinates (HETATM + ATOMS)
        10 Num. of CONECT
        11 Num. of SEQRES
    """
    master_widths = (6, 9, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5)
    master_parse = get_parser(master_widths)
    try:
        master_fields = master_parse(master_line)
    except UnboundLocalError:
        master_fields = None
        print(f"Warning! No MASTER line found in '{fn}'!")

    atoms = list()
    coords = list()
    fragments = {}

    for i, line in enumerate(atm_lines):
        # Charge is not considered right now ...
        fields = atm_parse(f"{line[:78]: <78}")
        res_name = fields[4].strip()
        chain = fields[5].strip()
        res_seq = int(fields[6])
        frag = f"{chain}_{res_seq}_{res_name}"

        xyz = fields[8:11]
        atom = fields[13].strip()
        if not atom.lower() in KNOWN_ATOMS:
            name = fields[2]
            atom = parse_atom_name(name)
        atoms.append(atom)
        coords.append(xyz)

        frag = fragments.setdefault(frag, list())
        frag.append(i)

    if master_fields is not None:
        # Verification using MASTER record
        num_coord = int(master_fields[9])
        assert len(atoms) == num_coord

    coords = np.array(coords, dtype=float) * ANG2BOHR
    return atoms, coords.flatten(), fragments


def geom_from_pdb(fn, **kwargs):
    atoms, coords, fragments = parse_pdb(fn)
    kwargs["fragments"] = fragments
    geom = Geometry(atoms, coords.flatten(), **kwargs)
    return geom


def get_conect_lines(atoms, coords):
    bonds = find_bonds(atoms, coords)
    # PDB indexing is 1-based
    bonds += 1
    # Bring in a suitable order for CONECT entries
    bonds.sort(axis=1)
    bonds = sorted(bonds, key=lambda row: row[0])
    conect = OrderedDict()
    for from_, to_ in bonds:
        conect.setdefault(from_, list()).append(to_)

    conect_lines = list()
    for from_, to_ in conect.items():
        to_iter = chunks(sorted(to_), 4)
        for to_chunk in to_iter:
            conect_lines.append([from_] + to_chunk)
    return conect_lines


def atoms_coords_to_pdb_str(atoms, coords, fragments=None, resname="", conect=True):
    coords3d = coords.reshape(-1, 3)
    coords3d_ang = coords3d / ANG2BOHR

    coord_fmt = "{: >8.3f}" * 3
    #            serial  name  altLoc      chainID     iCode
    #                                resName     resSeq
    hetatm_fmt = "HETATM{: >5d} {: >4}{: >1}{: >3}{: >1} {: >4d}{: >1}   "
    #      xyz
    hetatm_fmt += coord_fmt
    #      occupancy tempFactor            atom
    hetatm_fmt += "{: >6.2f}{: >6.2f}" + 10 * " " + "{: >2s}"
    ter_fmt = "TER   {: >5d}      {: >3s} {:1s}{: >4d}"

    if fragments is None:
        fragments = [
            range(len(atoms)),
        ]

    # Fixed for now
    altLoc = ""
    resName = resname
    chainID = ""
    iCode = ""
    occupancy = 1.0
    tempFactor = 0.0

    lines = list()
    serial = 1
    for resSeq, fragment in enumerate(fragments, 1):
        for id_ in fragment:
            name = atoms[id_]
            xyz = coords3d_ang[id_]
            line = hetatm_fmt.format(
                serial,
                name,
                altLoc,
                resName,
                chainID,
                resSeq,
                iCode,
                *xyz,
                occupancy,
                tempFactor,
                name,
            )
            lines.append(line)
            serial += 1
        lines.append(ter_fmt.format(serial, resName, chainID, resSeq))

    def fmt_conect_line(conect_line):
        return "CONECT" + ("{: >5d}" * len(conect_line)).format(*conect_line)

    if conect:
        conect_lines = get_conect_lines(atoms, coords)
    else:
        conect_lines = []
    conect_str = "\n".join([fmt_conect_line(cl) for cl in conect_lines])

    pdb_tpl = Template(
        textwrap.dedent(
            """\
    REMARK 1 Created by pysisyphus
    {%+ for line in lines -%}
    {{ line }}
    {%+ endfor -%}
    {{- conect_str }}
    END"""
        )
    )

    pdb_str = pdb_tpl.render(lines=lines, conect_str=conect_str)
    return pdb_str


def geom_to_pdb_str(geom, detect_fragments=False, **kwargs):
    fragments = None
    if detect_fragments:
        try:
            fragments = geom.internal.fragments
        except AttributeError:
            geom_ = geom.copy(coord_type="redund")
            fragments = geom_.internal.fragments
    pdb_str = atoms_coords_to_pdb_str(
        geom.atoms, geom.coords3d, fragments=fragments, **kwargs
    )
    return pdb_str
