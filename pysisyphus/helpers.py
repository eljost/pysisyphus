import os
from pathlib import Path

from pysisyphus.constants import ANG2BOHR
from pysisyphus.Geometry import Geometry
from qchelper.geometry import parse_xyz_file, parse_trj_file


def geom_from_xyz_file(xyz_fn):
    atoms, coords = parse_xyz_file(xyz_fn)
    coords *= ANG2BOHR
    geom = Geometry(atoms, coords.flatten())
    return geom


def geom_from_library(xyz_fn):
    this_dir = os.path.abspath(os.path.dirname(__file__))
    xyz_dir = Path(this_dir) / "../xyz_files/"
    atoms, coords = parse_xyz_file(xyz_dir / xyz_fn)
    coords *= ANG2BOHR
    geom = Geometry(atoms, coords.flatten())
    return geom


def geoms_from_trj(trj_fn):
    geoms = [Geometry(atoms, coords.flatten()*ANG2BOHR)
             for atoms, coords in parse_trj_file(trj_fn)
    ]
    return geoms


if __name__ == "__main__":
    print(load_geometry("hcn.xyz"))
    print(geoms_from_trj("cycle_040.trj"))
