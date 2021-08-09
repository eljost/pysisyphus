import numpy as np

from pysisyphus.constants import ANG2BOHR
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import file_or_str
from pysisyphus.xyzloader import split_xyz_str


@file_or_str(".xyz", ".trj")
def parse_xyz(xyz_str, with_comment=False):
    lines = iter(xyz_str.strip().split("\n"))
    atoms_coords = list()
    comments = list()
    for line in lines:
        line = line.strip()
        # Skip emtpy lines
        if not line:
            continue
        # A new xyz block starts with an integer
        try:
            expect = int(line)
        # We only skip empty lines, otherwise we leave the loop.
        except ValueError:
            break
        atoms = list()
        coords3d = np.zeros((expect, 3), dtype=float)
        # Advance iterator over xyz definition
        comment = next(lines)
        comments.append(comment if with_comment else None)
        for i in range(expect):
            atom, *xyz_ = next(lines).split()
            atoms.append(atom)
            coords3d[i] = xyz_
        assert i == expect - 1
        assert len(atoms) == expect
        coords3d *= ANG2BOHR
        atoms_coords.append((atoms, coords3d.flatten()))
    return atoms_coords, comments



def geoms_from_xyz(fn, **kwargs):
    atoms_coords, comments = parse_xyz(fn, with_comment=True)
    geoms = [
        Geometry(atoms, coords.flatten(), comment=comment, **kwargs)
        for (atoms, coords), comment in zip(atoms_coords, comments)
    ]
    return geoms


def geom_from_xyz(fn, **kwargs):
    return geoms_from_xyz(fn, **kwargs)[0]


def geoms_from_inline_xyz(inline_xyz, **kwargs):
    atoms_coords = split_xyz_str(inline_xyz)
    # We excpect the coordinates to be given in Angstrom
    geoms = [
        Geometry(atoms, coords * ANG2BOHR, **kwargs) for atoms, coords in atoms_coords
    ]
    return geoms
