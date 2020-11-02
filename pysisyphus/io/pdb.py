import struct
import textwrap

from jinja2 import Template
import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.constants import ANG2BOHR


AMINOS = "ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU" \
         "MET ASN PRO GLN ARG SER THR VAL TRP TYR".split()


def get_parser(widths):
    fmt = " ".join('{}{}'.format(width, "s") for width in widths)
    fieldstruct = struct.Struct(fmt)
    parse = lambda line: tuple(s.decode() for s in fieldstruct.unpack(line.encode()))
    return parse


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
        fields = atm_parse(line[:78])
        res_name = fields[4].strip()
        chain = fields[5].strip()
        res_seq = int(fields[6])
        frag = f"{chain}_{res_seq}_{res_name}"

        xyz = fields[8:11]
        atom = fields[13].strip()
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


def atoms_coords_to_pdb_str(atoms, coords, fragments=None):
    coords3d = coords.reshape(-1, 3)

    coord_fmt = "{: >8.3f}"*3
    #            serial  name  altLoc      chainID     iCode
    #                                resName     resSeq           
    fmt = "HETATM{: >5d} {: >4}{: >1}{: >3}{: >1} {: >4d}{: >1}   "
    #      xyz
    fmt += coord_fmt
    #      occupancy tempFactor            atom
    fmt += "{: >6.2f}{: >6.2f}" + 10*" " + "{: >2s}"

    if fragments is None:
        fragments = [range(len(atoms)), ]

    # Fixed for now
    altLoc = ""
    resName = ""
    chainID = ""
    iCode = ""
    occupancy = 1.0
    tempFactor = 0.0

    lines = list()
    serial = 0
    for resSeq, fragment in enumerate(fragments):
        for id_ in fragment:
            name = atoms[id_]
            xyz = coords3d[id_]
            line = fmt.format(serial, name, altLoc, resName, chainID, resSeq,
                              iCode, *xyz, occupancy, tempFactor, name)
            lines.append(line)
            serial += 1

    pdb_tpl = Template(textwrap.dedent(
    """\
    REMARK 1 Created by pysisyphus
    {%+ for line in lines -%}
    {{ line }}
    {%+ endfor -%} 
    END"""))

    pdb_str = pdb_tpl.render(lines=lines)
    return pdb_str


def geom_to_pdb_str(geom, detect_fragments=False):
    # Convert to Å
    coords = geom.cart_coords / ANG2BOHR
    if len(geom.fragments) > 0:
        raise Exception("Not yet implemented!")

    fragments = None
    if detect_fragments:
        try:
            fragments = geom.internal.fragments
        except AttributeError:
            geom_ = geom.copy(coord_type="redund")
            fragments = geom_.internal.fragments
    pdb_str = atoms_coords_to_pdb_str(geom.atoms, coords, fragments=fragments)
    return pdb_str
