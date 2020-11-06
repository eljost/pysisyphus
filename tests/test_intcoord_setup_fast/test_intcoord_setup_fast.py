from pysisyphus.helpers import geom_loader
from pysisyphus.intcoords.setup_fast import find_bonds, find_bonds_bends


def test_find_bonds_fast():
    geom = geom_loader("lib:split.image_021.xyz")

    bonds = find_bonds(geom)
    assert len(bonds) == 31


def test_find_bonds_bends_fast():
    geom = geom_loader("lib:split.image_021.xyz")

    bonds, bends = find_bonds_bends(geom)
    assert len(bonds) == 31
    assert len(bends) == 51
