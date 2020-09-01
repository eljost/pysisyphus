import itertools as it

import numpy as np
import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.intcoords.helpers import form_coordinate_union
from pysisyphus.intcoords.RedundantCoords import RedundantCoords
from pysisyphus.run import run_from_dict
from pysisyphus.testing import using


def test_coordinate_union():
    coord_type = "redund_v2"
    geoms = geom_loader("rx_001.trj", coord_type=coord_type)
    ts = geom_loader("ts_001.xyz", coord_type=coord_type)
    print(geoms)
    pics = list()
    for geom in (*geoms, ts):
        int_ = geom.internal
        pic = int_.prim_indices_set
        print(len(int_._primitives))
        pics.append(pic)

    for i, p1 in enumerate(pics):
        for j, p2 in enumerate(pics[i+1:]):
            union = p1 | p2
            intersection = p1 & p2
            ul = len(union)
            il = len(intersection)
            print(f"({i},{j}): p1={len(p1)}, p2={len(p2)}, union={ul}, intersection={il}")

    geom1, geom2 = geoms
    prim_indices = form_coordinate_union(geom1, geom2)
    # s = sum([len(p) for p in prim_indices])
    # print(s)
    ts_ = RedundantCoords(ts.atoms, ts.cart_coords, prim_indices=prim_indices)
    ts__ = geom_loader("ts_001.xyz", coord_type="redund_v2",
                       coord_kwargs={"prim_indices": prim_indices,})
    import pdb; pdb.set_trace()


@using("pyscf")
def test_run_geom_section_union():
    ref_energy = -160.7433945831947

    def run_assert(run_dict, prim_len):
        res = run_from_dict(run_dict)
        geom = res.calced_geoms[0]
        assert geom.energy == pytest.approx(ref_energy)
        int_ = geom.internal
        assert len(int_.prim_indices_set) == prim_len

        # H = geom.hessian
        # w, v = np.linalg.eigh(H)
        # first = 5
        # print("eigvals", w[:first])

    def get_run_dict():
        return {
            "calc": {
                "type": "pyscf",
                "basis": "sto3g",
            },
            "geom": {
                "type": "redund_v2",
                "fn": "ts_001.xyz",
            },
        }

    run_assert(get_run_dict(), prim_len=56)

    run_dict = get_run_dict()
    run_dict["geom"]["union"] = "rx_001.trj"
    run_assert(run_dict, prim_len=66)
