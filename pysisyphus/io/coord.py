import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.helpers_pure import file_or_str


@file_or_str("coord", exact=True)
def geom_from_turbomole_coord(text, **geom_kwargs):
    text = text.strip()
    assert text.startswith("$coord")
    lines = text.split("\n")
    first_line = lines.pop(0)
    assert "$coord" in first_line
    xyzs = list()
    atoms = list()
    for line in lines:
        line = line.strip()
        # Only parse until the first $ after $coord is encountered, e.g.,
        # when the $redundant section starts
        if line.startswith("$"):
            break
        *xyz, atom = line.split()
        xyzs.append(xyz)
        atoms.append(atom)
    coords = np.array(xyzs, dtype=float).flatten()
    geom = Geometry(atoms, coords, **geom_kwargs)
    return geom


@file_or_str("coord", exact=True)
def geoms_from_turbomole_coord(text, **geom_kwargs):
    geom = geom_from_turbomole_coord(text, **geom_kwargs)
    return (geom,)
