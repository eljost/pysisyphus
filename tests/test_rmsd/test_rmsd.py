import numpy as np
from numpy.random import MT19937, RandomState, SeedSequence
from pytest import approx

from pysisyphus.linalg import rmsd_grad


def test_rmsd_grad():
    rs = RandomState(MT19937(SeedSequence(26062020)))
    n = 3
    coords3d = np.eye(3)
    ref_coords3d = coords3d + 0.2 * rs.random((n, 3))

    rmsd, grad = rmsd_grad(coords3d, ref_coords3d)

    assert rmsd == approx(0.05086180)
    ref_grad = np.array(
        [
            [0.29974081, -0.07879816, -0.22094265],
            [-0.07879816, -0.15494155, 0.23373971],
            [-0.22094265, 0.23373971, -0.01279706],
        ]
    )
    np.testing.assert_allclose(grad, ref_grad, atol=1e-8)
