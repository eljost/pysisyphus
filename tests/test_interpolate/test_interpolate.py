#!/usr/bin/env python3

from pysisyphus.helpers import geom_from_library, geom_from_xyz_file
from pysisyphus.interpolate.LST import LST
from pysisyphus.interpolate.IDPP import IDPP


def test_lst():
    initial = geom_from_library("dipeptide_init.xyz")
    final = geom_from_library("dipeptide_fin.xyz")

    geoms = (initial, final)
    lst = LST(geoms, 50, align=True)
    geoms = lst.interpolate_all()
    lst.all_geoms_to_trj("lst_opt.trj")


def test_idpp():
    # initial = geom_from_library("dipeptide_init.xyz")
    # final = geom_from_library("dipeptide_fin.xyz")

    initial = geom_from_xyz_file("09_htransfer_product.xyz")
    final = geom_from_xyz_file("10_po_diss_product_xtbopt.xyz")

    geoms = (initial, final)
    idpp = IDPP(geoms, 18, align=True)
    geoms = idpp.interpolate_all()
    idpp.all_geoms_to_trj("idpp_opt.trj")


if __name__ == "__main__":
    # test_lst()
    test_idpp()
