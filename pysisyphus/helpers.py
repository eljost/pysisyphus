import os
from pathlib import Path

from pysisyphus.Geometry import Geometry
from qchelper.geometry import parse_xyz_file


def load_geometry(xyz_fn):
    this_dir = os.path.abspath(os.path.dirname(__file__))
    xyz_dir = Path(this_dir) / "../xyz_files/"
    atoms, coords = parse_xyz_file(xyz_dir / xyz_fn)
    geom = Geometry(atoms, coords.flatten())
    return geom


if __name__ == "__main__":
    load_geometry("hcn.xyz")
