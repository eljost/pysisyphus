import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.io.pdb import geom_to_pdb_str


@pytest.mark.parametrize(
    "pdb_fn, fragment_num", [
        ("lib:pdbs/1gcn.pdb", 29),
        ("lib:pdbs/1bl8.pdb", 388+4),
    ]
)
def test_fragment_num(pdb_fn, fragment_num):
    geom = geom_loader(pdb_fn)
    # geom.jmol()

    assert len(geom.fragments) == fragment_num


def test_get_fragments():
    full_geom = geom_loader("lib:pdbs/1bl8.pdb")
    geom = full_geom.get_fragments("75_THR")
    # geom.jmol()

    assert len(geom.fragments) == 4


def test_pdb_write(this_dir):
    geom = geom_loader("lib:h2o.xyz")
    pdb_str = geom_to_pdb_str(geom)

    # with open("h2o_ref.pdb", "w") as handle:
        # handle.write(pdb_str)

    # Reference pdb
    with open(this_dir / "h2o_ref.pdb") as handle:
        ref = handle.read()

    assert pdb_str == ref


def test_geom_to_pdb(this_dir):
    geom = geom_loader(this_dir / "five_chloroforms_xtbopt.xyz")
    pdb_str = geom_to_pdb_str(geom, detect_fragments=True)

    with open("five_chloroforms_ref.pdb", "w") as handle:
        handle.write(pdb_str)

    # Reference pdb
    with open(this_dir / "five_chloroforms_ref.pdb") as handle:
        ref = handle.read()

    assert pdb_str == ref
