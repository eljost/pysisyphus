from pysisyphus.helpers import geom_loader
from pysisyphus.intcoords.setup_fast import (
    find_bonds_for_geom,
    find_bonds_bends,
    get_bond_vec_getter,
)


def test_find_bonds_fast():
    geom = geom_loader("lib:split.image_021.xyz")

    bonds = find_bonds_for_geom(geom)
    assert len(bonds) == 31


def test_find_bonds_bends_fast():
    geom = geom_loader("lib:split.image_021.xyz")

    bonds, bends = find_bonds_bends(geom)
    assert len(bonds) == 31
    assert len(bends) == 51


def test_bond_vec_getter():
    geom = geom_loader("lib:acetaldehyd_oniom.xyz")
    bond_vec_getter = get_bond_vec_getter(geom.atoms, geom.covalent_radii, (0, ))
    bond_vecs = bond_vec_getter(geom.coords3d)[0]
    assert len(bond_vecs) == 4
