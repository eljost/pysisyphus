import numpy as np

from pysisyphus.constants import ANG2BOHR
from pysisyphus.Geometry import Geometry
from pysisyphus.intcoords import RedundantCoords
from pysisyphus.intcoords.setup import setup_redundant
from pysisyphus.helpers_pure import get_cubic_crystal


def get_chain(num):
    mono_atoms = ["C", "C", "F", "H", "O", "H", "H"]
    CC = 1.5
    CH = 1.09
    CF = 1.39
    CO = 1.35
    mono_coords3d = np.array(
        (
            (0.0, 0.0, 0.0),  # C
            (CC, 0.0, 0.0),  # C
            (0.0, CF, 0.0),  # F
            (0.0, -CH, 0.0),  # H
            (CC, CO, 1.0),  # O
            (CC / 2, 1.5 * CO, 1.0),  # H-O
            (CC, -CH, 0.0),  # H
        )
    )
    mono_coords3d *= ANG2BOHR
    atoms = mono_atoms * num
    coords3d = np.tile(mono_coords3d, (num, 1))
    translation = np.arange(num) * 2 * CC * ANG2BOHR
    flip = np.ones(num)
    flip[1::2] = -1
    # Repeat along x-axis
    coords3d.reshape(num, -1, 3)[:, :, 0] += translation[:, None]
    # Mirror on y-axis
    # coords3d.reshape(num, -1, 3)[:, :, 1] *= flip[:, None]
    geom = Geometry(atoms, coords3d)
    return geom


# Or repeat H2O2 molecules in a chain
# Profile
#   python -m cProfile -o prof.pyprof `pyenv which pysistrj` num3.xyz --internals
# KCachegrind
#   pyprof2calltree -i prof.pyprof -k 
# Reuse KDTree?


def test_setup():
    for num in range(1, 10):
        geom = get_chain(num)
        # geom.jmol()
        with open(f"num{num}.xyz", "w") as handle:
            handle.write(geom.as_xyz())
        internal = RedundantCoords(geom.atoms, geom.coords3d)
        tps = len(internal.typed_prims)
        print(f"{num=}, {tps=}")


def test_cubic_crystal():
    atoms, coords3d = get_cubic_crystal(8.3, n=1)
    coord_info = setup_redundant(atoms, coords3d)
    assert len(coord_info.bonds) == 0
    assert len(coord_info.interfrag_bonds) == 54
