import struct

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
    master_fields = master_parse(master_line)

    atoms = list()
    coords = list()
    fragments = {}

    for i, line in enumerate(atm_lines):
        fields = atm_parse(line)
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

    # Verification using MASTER record
    num_coord = int(master_fields[9])
    assert len(atoms) == num_coord

    coords = np.array(coords, dtype=float) * ANG2BOHR
    return atoms, coords.flatten(), fragments


def geom_from_pdb(fn):
    atoms, coords, fragments = parse_pdb(fn)
    geom = Geometry(atoms, coords.flatten(), fragments=fragments)
    return geom
