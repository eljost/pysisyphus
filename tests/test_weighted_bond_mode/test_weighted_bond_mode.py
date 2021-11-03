import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.intcoords.helpers import get_weighted_bond_mode_getter


def test_get_weighted_bond_mode():
    geoms = geom_loader("lib:test_bond_mode.trj")
    #         break C-O    break N-C     form C-N
    target = [(9, 53, -1), (11, 54, -1), (9, 11, 1)]
    get_weighted_bond_mode = get_weighted_bond_mode_getter(target)
    masks = (
        (1, 1, 1),  # All targets are still valid at the first step
        (0, 1, 1),  # 9-53 is broken in all subsequent steps
        (0, 1, 1),
        (0, 1, 1),
        (0, 0, 0),  # 9-11 only broken in the last step
    )
    for mask, geom in zip(masks, geoms):
        bonds = get_weighted_bond_mode(geom.atoms, geom.coords3d)
        ref_bonds = [b for b, m in zip(target, mask) if m]
        assert bonds == ref_bonds
        print(mask)


@pytest.mark.skip
def test_get_frac_weighted_bond_mode(this_dir):
    geom = geom_loader(this_dir / "01_input.xyz")
    #         break C-O    break N-C     form C-N
    # target = [(9, 53, -1), (11, 54, -1), (9, 11, 1)]
    target = ((53, 9, -1), (11, 54, -1), (9, 11, 1), (8, 34, 1))
    print()
    get_weighted_bond_mode = get_weighted_bond_mode_getter(target, fractional=True)
    bonds = get_weighted_bond_mode(geom.atoms, geom.coords3d)
    print(bonds)
