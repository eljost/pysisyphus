import numpy as np

from pysisyphus.constants import ANG2BOHR
from pysisyphus.helpers_pure import file_or_str


@file_or_str(".xyz", ".trj")
def parse_xyz_v2(xyz_str, with_comment=False):
    lines = iter(xyz_str.strip().split("\n"))
    # geoms = list()
    atoms_coords = list()
    for l in lines:
        l = l.strip()
        # Skip emtpy lines
        if not l:
            continue
        # A new xyz block starts with an integer
        try:
            expect = int(l)
        except ValueError:
            continue
        atoms = list()
        coords3d = np.zeros((expect, 3), dtype=float)
        # Advance iterator over xyz definition
        comment = next(lines)
        comment = comment if with_comment else None
        for i in range(expect):
            atom, *xyz_ = next(lines).split()
            atoms.append(atom)
            coords3d[i] = xyz_
        assert i == expect - 1
        assert len(atoms) == expect
        coords3d *= ANG2BOHR
        atoms_coords.append((atoms, coords3d.flatten()))
    return atoms_coords
