import pytest

from pysisyphus.io.pubchem import cid_from_name, sdf_from_cid, geom_from_pubchem_name


@pytest.mark.skip
def test_pubchem_cid_for_name():
    name = "sodiumborohydride"
    cid = cid_from_name(name)
    assert cid == 4311764


@pytest.mark.skip
@pytest.mark.parametrize(
    "name", [
        "benzene",
        "methane",
    ]
)
def test_geom_from_pubchem_name(name):
    geom = geom_from_pubchem_name(name)
    assert geom
    # geom.jmol()
