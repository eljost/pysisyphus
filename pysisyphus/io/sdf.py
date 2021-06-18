import numpy as np

from pysisyphus.constants import BOHR2ANG
from pysisyphus.Geometry import Geometry


def parse_coord_line(line):
    x, y, z, atom, *_ = line.strip().split()
    return (x, y, z), atom


def parse_sdf(text):
    lines = text.split("\n")
    # title, program, comment = lines[:3]
    count = lines[3]
    atoms, *_ = count.split()
    atoms = int(atoms)

    coord_lines = lines[4:4+int(atoms)]
    coords, atoms = zip(*[parse_coord_line(cl) for cl in coord_lines])
    coords = np.array(coords, dtype=float)
    coords /= BOHR2ANG
    return atoms, coords


def geom_from_sdf(text, **kwargs):
    atoms, coords = parse_sdf(text)
    geom = Geometry(atoms, coords, **kwargs)
    return geom
