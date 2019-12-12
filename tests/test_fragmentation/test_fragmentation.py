#!/usr/bin/env python3

import numpy as np

from pysisyphus.helpers import geom_from_xyz_file
from pysisyphus.constants import BOHR2ANG
from pysisyphus.xyzloader import make_trj_str


def fragment(geom, prim_coord, step=6, step_num=50):
    assert geom.coord_type != "cart"

    index_ = geom.internal.get_index_of_prim_coord(prim_coord)
    prim_val = geom.coords[index_]
    print(prim_coord, index_, prim_val)

    step_ = step / step_num
    base_step = np.zeros_like(geom.coords)
    base_step[index_] = step / step_num

    cart_coords = list()
    for i in range(step_num):
        print(f"{i:02d} ... ", end="")
        new_coords = geom.coords + base_step
        geom.coords = new_coords
        print(not geom.internal.backtransform_failed)
        cart_coords.append(geom.cart_coords.copy())
    return cart_coords


def test_fragmentation():
    fn = "01_dihedral_scan.009.xyz"
    geom = geom_from_xyz_file(fn, coord_type="redund")

    prim_coord = (19, 20)
    bbis = geom.internal.bare_bond_indices.tolist()
    _ = bbis.copy()
    bbis.remove(list(sorted(prim_coord)))
    import pdb; pdb.set_trace()
    # bbis.delete(

    frag_kwargs = {
        "step": 6,
        "step_num": 10,
    }
    cc = fragment(geom, prim_coord, **frag_kwargs)
    cc = np.reshape(cc, (-1, len(geom.atoms), 3)) * BOHR2ANG
    trj_str = make_trj_str(geom.atoms, cc)
    with open("fragmented.xyz", "w") as handle:
        handle.write(trj_str)
    # Easier fragment detection: just delete the supplied prim_coord ...


if __name__ == "__main__":
    test_fragmentation()
