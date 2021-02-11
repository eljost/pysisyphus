import pytest

from pysisyphus.io import geoms_from_molden, geoms_from_xyz


def test_molden(this_dir):
    molden_fn = this_dir / "geometries.molden"
    xyz_fn = this_dir / "ref.xyz"

    molden_geoms = geoms_from_molden(molden_fn)
    ref_geoms = geoms_from_xyz(xyz_fn)

    for mg, rg in zip(molden_geoms, ref_geoms):
        assert mg.rmsd(rg) == pytest.approx(0.0)

    assert len(molden_geoms) == 5
    assert len(mg.atoms) == 28
