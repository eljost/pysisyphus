import pytest

from pysisyphus.io.pubchem import cid_for_name, sdf_for_cid, geom_for_pubchem_name


def test_pubchem_cid_for_name():
    name = "sodiumborohydride"
    cid = cid_for_name(name)
    assert cid == 4311764


@pytest.mark.parametrize(
    "name", [
        "benzene",
        "methane",
    ]
)
def test_geom_from_pubchem_name(name):
    geom = geom_for_pubchem_name(name)
    assert geom
    # geom.jmol()
