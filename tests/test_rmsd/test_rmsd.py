import numpy as np
from numpy.random import MT19937, RandomState, SeedSequence
from pytest import approx

from pysisyphus.finite_diffs import finite_difference_gradient
from pysisyphus.linalg import fd_rmsd_hessian, rmsd_grad


def test_rmsd_grad():
    rs = RandomState(MT19937(SeedSequence(26062020)))
    n = 3
    coords3d = np.eye(3)
    ref_coords3d = coords3d + 0.2 * rs.random((n, 3))

    rmsd, grad = rmsd_grad(coords3d, ref_coords3d)

    def scalar_func(coords):
        coords3d = coords.reshape(-1, 3)
        return rmsd_grad(coords3d, ref_coords3d)[0]

    fd_grad = finite_difference_gradient(
        coords3d.flatten(), scalar_func, step_size=1e-3
    )
    fd_grad3d = fd_grad.reshape(-1, 3)
    # Test gradient against finite differences
    np.testing.assert_allclose(grad, fd_grad3d, atol=5e-6)

    # Assert that everything is at a certain value
    assert rmsd == approx(0.05086180)
    ref_grad = np.array(
        [
            [0.29974081, -0.07879816, -0.22094265],
            [-0.07879816, -0.15494155, 0.23373971],
            [-0.22094265, 0.23373971, -0.01279706],
        ]
    )
    np.testing.assert_allclose(grad, ref_grad, atol=1e-8)


def test_rmsd_hessian():
    rs = RandomState(MT19937(SeedSequence(26062020)))
    n = 3
    coords3d = np.eye(3)
    ref_coords3d = coords3d + 0.2 * rs.random((n, 3))
    rmsd, hessian = fd_rmsd_hessian(coords3d, ref_coords3d)
    assert rmsd == approx(0.05086180)
    np.testing.assert_allclose(hessian, hessian.T)
    w, _ = np.linalg.eigh(hessian)
    assert w[0] == approx(-3.36876688e-01)
