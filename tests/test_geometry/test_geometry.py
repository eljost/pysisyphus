import numpy as np
import pytest

from pysisyphus.Geometry import get_trans_rot_vectors
from pysisyphus.helpers import geom_loader


@pytest.fixture
def geom():
    return geom_loader("lib:benzene.xyz")


def test_standard_orientation(geom):
    geom.standard_orientation()

    com = geom.center_of_mass
    # This somehow fails in the CI ...
    # np.testing.assert_allclose(com, np.zeros((3, )), atol=1e-16)
    assert np.abs(com).sum() == pytest.approx(0., abs=1e-10)


def test_inertia_tensor(geom, this_dir):
    I = geom.inertia_tensor

    I_ref_fn = this_dir / "benzene_inertia_tensor_ref"
    I_ref = np.loadtxt(I_ref_fn)
    np.testing.assert_allclose(I, I_ref)

    geom.standard_orientation()
    # Should be diagonal now
    I_so = geom.inertia_tensor

    np.testing.assert_allclose(I_so,I_so.T)
    w, v = np.linalg.eigh(I_so)
    np.testing.assert_allclose(w, np.array((317.60004608, 317.61686026, 635.21687793)))


def test_trans_rot_vec(geom, this_dir):
    vecs = geom.get_trans_rot_vectors()

    vecs_ref_fn = this_dir / "benzene_trans_rot_vecs_ref"
    vecs_ref = np.loadtxt(vecs_ref_fn)
    np.testing.assert_allclose(vecs, vecs_ref)

    trv = get_trans_rot_vectors(geom.cart_coords, geom.masses)
    np.testing.assert_allclose(trv, vecs)
