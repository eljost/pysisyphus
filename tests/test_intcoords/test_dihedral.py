import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.elem_data import COVALENT_RADII as CR, MASS_DICT as MD
from pysisyphus.trj import print_internals


def test_dihedral_definition():
    atoms = "Q Q Q Q".split()
    CR["q"] = 0.5
    MD["q"] = 1
    coords = np.array((
                (0, 0, 0),
                (1, 0, 0),
                (1, 1, 0),
                (0, 1, 0),
    ), dtype=float)

    geom = Geometry(atoms, coords.flatten())
    # print(geom)
    # print(geom.coords)

    geom2 = geom.copy(coord_type="redund")
    # print(geom2)
    # print(geom2.coords)

    print_internals((geom, ))
    int_ = geom2.internal
    prim_inds = set([tuple(pi.inds) for pi in int_._prim_internals])

    ref = {
        # 4 Bonds
        (0, 1),
        (1, 2),
        (2, 3),
        (0, 3),
        # 4 Bends
        (0, 1, 2),
        (1, 2, 3),
        (0, 3, 2),
        (1, 0, 3),
        # 4 Dihedrals
        (1, 0, 3, 2),
        (0, 3, 2, 1),
        (2, 1, 0, 3),
        (0, 1, 2, 3),
    }

    assert (prim_inds - ref) == set()
