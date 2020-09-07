import pytest
from rmsd import kabsch_rmsd

from pysisyphus.constants import ANG2BOHR as AB
from pysisyphus.helpers import geom_loader
from pysisyphus.io.zmat import geom_from_zmat, zmat_from_fn, ZLine


def assert_geom(ref_fn, zmat_fn, atol=2.5e-5):
    zmat = zmat_from_fn(zmat_fn)
    geom = geom_from_zmat(zmat)

    ref = geom_loader(ref_fn)
    rmsd = kabsch_rmsd(geom.coords3d, ref.coords3d, translate=True)
    print(f"RMSD: {rmsd:.6f}")
    assert rmsd == pytest.approx(0., abs=atol)


def test_glycine_resorted(this_dir):
    assert_geom(this_dir / "glycine_resorted.xyz",
                this_dir / "glycine_resorted.zmat")


def test_geom_from_zmat(this_dir):
    zmat = [
        ZLine("C"),
        ZLine("C", 0, 1.510486*AB),
        ZLine("N", 0, 1.459785*AB, 1, 112.257683),
        ZLine("O", 1, 1.220389*AB, 0, 118.653885, 2, -179.541229),
        ZLine("O", 1, 1.353023*AB, 0, 122.591621, 2, 0.825231),
    ]

    geom = geom_from_zmat(zmat)

    reference = geom_loader(this_dir / "glycine_noh.xyz")
    rmsd = kabsch_rmsd(geom.coords3d, reference.coords3d, translate=True)
    assert rmsd == pytest.approx(0., abs=1e-6)
