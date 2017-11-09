import os
from pathlib import Path

from pysisyphus.Geometry import Geometry
from qchelper.geometry import parse_xyz_file, parse_trj_file


def geom_from_xyz_file(xyz_fn):
    atoms, coords = parse_xyz_file(xyz_fn)
    geom = Geometry(atoms, coords.flatten())
    geom.set_calculator(CALC_DICT[args.calc]())
    return geom


def geom_from_library(xyz_fn):
    this_dir = os.path.abspath(os.path.dirname(__file__))
    xyz_dir = Path(this_dir) / "../xyz_files/"
    atoms, coords = parse_xyz_file(xyz_dir / xyz_fn)
    geom = Geometry(atoms, coords.flatten())
    return geom


def geoms_from_trj(trj_fn):
    print(trj_fn)
    geoms = [Geometry(atoms, coords.flatten())
             for atoms, coords in parse_trj_file(trj_fn)
    ]
    return geoms


if __name__ == "__main__":
    #load_geometry("hcn.xyz")
    geoms_from_trj("cycle_040.trj")
