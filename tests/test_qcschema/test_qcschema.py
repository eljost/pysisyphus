import json

import pytest

from pysisyphus.io import geom_from_qcschema


@pytest.fixture
def qcschema_str():
    return """{
        "schema_name": "qc_schema_input",
        "schema_version": 1,
        "molecule": {
            "geometry": [0.0, 0.0, -0.1294, 0.0, -1.4941, 1.0274, 0.0, 1.4941, 1.0274],
            "symbols": ["O", "H", "H"]
        }
    }"""


@pytest.fixture
def qcschema_dict(qcschema_str):
    return json.loads(qcschema_str)


def test_geom_from_qcschema_str(qcschema_str, qcschema_dict):
    geom1 = geom_from_qcschema(qcschema_str)
    geom2 = geom_from_qcschema(qcschema_dict)
    assert geom1 == geom2


def test_geom_from_qcschema_kwargs(qcschema_str):
    geom_kwargs = {
            "coord_type": "redund",
    }
    geom = geom_from_qcschema(qcschema_str, **geom_kwargs)
    assert geom.coord_type == "redund"
    assert len(geom.coords) == 3
