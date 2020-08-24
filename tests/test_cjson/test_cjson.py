import numpy as np

from pysisyphus.io import geom_from_cjson


def test_cjson(this_dir):
    fn = this_dir / "test.cjson"
    geom = geom_from_cjson(fn)
    np.testing.assert_allclose(geom.centroid, (-3.92621311, 0.95627228, -0.03928268))
    assert len(geom.atoms) == 8
    # geom.jmol()
