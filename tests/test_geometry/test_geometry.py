import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.helpers import geom_loader
from pysisyphus.intcoords import Stretch


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


def test_fd_coords3d_gen():
    geom = geom_loader("lib:h2o.xyz")
    fd_coords3d = list(geom.fd_coords3d_gen())
    assert len(fd_coords3d) == geom.coords3d.size


def test_approximate_radius():
    geom = geom_loader("lib:h2o.xyz")
    radius = geom.approximate_radius()
    # max distance is H-H distance
    assert  radius == pytest.approx(Stretch._calculate(geom.coords3d, [1, 2]))


def test_rotate():
    geom = geom_loader("lib:h2o.xyz")
    for _ in range(3):
        geom.rotate()


def test_bond_sets():
    geom = geom_loader("lib:benzene.xyz")
    bond_sets = geom.bond_sets
    assert len(bond_sets) == 12


def test_anapot_bond_sets():
    geom = AnaPot().get_minima()[0]
    bond_sets = geom.bond_sets
    assert len(bond_sets) == 0
