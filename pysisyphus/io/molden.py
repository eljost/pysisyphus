import re

from pysisyphus.constants import ANG2BOHR
from pysisyphus.Geometry import Geometry
from pysisyphus.io.xyz import parse_xyz
from pysisyphus.helpers_pure import file_or_str


@file_or_str(".molden", ".input")
def parse_molden(molden_str):
    # molden_str = molden_str.strip()

    _, geometries = re.split("\[GEOMETRIES\] \(XYZ\)", molden_str)
    atoms_coords, comments = parse_xyz(geometries, with_comment=True)

    return atoms_coords, comments


def geoms_from_molden(fn, **kwargs):
    atoms_coords, comments = parse_molden(fn)
    geoms = [
        Geometry(atoms, coords.flatten(), comment=comment, **kwargs)
        for (atoms, coords), comment in zip(atoms_coords, comments)
    ]
    return geoms
