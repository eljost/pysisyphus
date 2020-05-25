import pytest

from pysisyphus.helpers import geom_loader


@pytest.mark.parametrize(
    "pdb_fn, fragment_num", [
        ("lib:pdbs/1gcn.pdb", 29),
        ("lib:pdbs/1bl8.pdb", 388+4),
    ]
)
def test_fragment_num(pdb_fn, fragment_num):
    geom = geom_loader(pdb_fn)
    # geom.jmol()

    assert len(geom.fragments) == fragment_num#29


def test_get_fragments():
    full_geom = geom_loader("lib:pdbs/1bl8.pdb")
    geom = full_geom.get_fragments("75_THR")
    # geom.jmol()

    assert len(geom.fragments) == 4
