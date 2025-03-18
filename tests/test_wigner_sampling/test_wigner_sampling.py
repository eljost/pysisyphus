import numpy as np
from scipy import stats

from pysisyphus.dynamics import wigner


def test_h2_wigner_sampling(this_dir):
    h5_fn = this_dir / "h2_hfdef2svp_final_hessian.h5"
    seed = 20182020
    sampler, _ = wigner.get_wigner_sampler(str(h5_fn), temperature=300, seed=seed)
    n = 500
    displ_coords3d = np.zeros((n, 2, 3))
    for i in range(n):
        displ_coords3d[i], _ = sampler()
    dists = np.linalg.norm(displ_coords3d[:, 0] - displ_coords3d[:, 1], axis=1)
    descr = stats.describe(dists)
    min_, max_ = descr.minmax
    mean = descr.mean
    var = descr.variance
    np.testing.assert_allclose(
        (min_, max_, mean, var), (1.01207, 1.88892, 1.41359, 0.02420), atol=1e-5
    )
