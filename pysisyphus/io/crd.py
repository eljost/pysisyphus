import re

import numpy as np

from pysisyphus.constants import ANG2BOHR
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import file_or_str


@file_or_str(".crd")
def parse_crd(crd):
    lines = [l.strip() for l in crd.strip().split("\n") if not l.startswith("*")]
    expect_line = lines.pop(0)
    expect_num = int(expect_line.split()[0])
    coords3d = np.zeros((expect_num, 3), dtype=float)
    atoms = list()
    re_ = re.compile(r"\d+")
    for i, line in enumerate(lines):
        _, _, res, atom, x, y, z, *_ = line.split()
        coords3d[i] = (x, y, z)
        atom = re_.sub("", atom)  # Delte number
        atoms.append(atom)
    assert len(atoms) == expect_num
    coords3d *= ANG2BOHR
    return atoms, coords3d.flatten()


def geom_from_crd(fn, **kwargs):
    atoms, coords = parse_crd(fn)
    return Geometry(atoms, coords, **kwargs)
