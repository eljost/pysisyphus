from pysisyphus.helpers import geom_loader
from pysisyphus.intcoords.setup import get_dihedral_inds
from pysisyphus.intcoords.setup_fast import (
    find_bonds_for_geom,
    find_bonds_bends,
    find_dihedrals,
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
    bond_vec_getter = get_bond_vec_getter(geom.atoms, geom.covalent_radii, (0,))
    bond_vecs = bond_vec_getter(geom.coords3d)[0]
    assert len(bond_vecs) == 4


def test_find_dihedrals():
    geom = geom_loader("lib:hydrogen_bond_fragments_test.xyz")
    max_deg = 175
    bonds, bends = find_bonds_bends(geom)
    ref_dihedrals, _ = get_dihedral_inds(geom.coords3d, bonds, bends, max_deg=max_deg)
    dihedrals = find_dihedrals(geom.coords3d, bonds, bends, max_deg=max_deg)
    # print(f" Ref dihedrals: {len(ref_dihedrals)}")
    # print(f"Fast dihedrals: {len(dihedrals)}")
    assert len(ref_dihedrals) == 24
    assert len(dihedrals) == 32
