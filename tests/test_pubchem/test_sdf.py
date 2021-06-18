import pytest

from pysisyphus.io.sdf import parse_sdf, geom_from_sdf


@pytest.fixture
def sdf_text(this_dir):
    fn = this_dir / "4311764.sdf"
    with open(fn) as handle:
        text = handle.read()
    return text


def test_parse_sdf(sdf_text):
    atoms, coords = parse_sdf(sdf_text)
    assert len(atoms) == 6
