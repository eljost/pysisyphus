from pathlib import Path
import tempfile

import numpy as np

from pysisyphus.io import geom_from_hessian
from pysisyphus.marcus_dim import fit_marcus_dim, get_batch_calc_func_from_npz


def test_marcus_dim_fit_replay(this_dir):
    """Replay fit of Marcus dimension from previous data stored in .npz file."""

    npz_fn = this_dir / "pdnb_xtb_marcus_dim_fit.npz"
    h5_fn = this_dir / "pdnb_xtb_final_hessian.h5"

    batch_calc_func, ncalcs = get_batch_calc_func_from_npz(npz_fn)
    geom = geom_from_hessian(h5_fn)
    with tempfile.TemporaryDirectory() as tmp_dir:
        out_dir = Path(tmp_dir)
        result = fit_marcus_dim(
            geom,
            batch_calc_func=batch_calc_func,
            batch_kwargs=dict(
                rms_thresh=0.005,
                batch_size=ncalcs,
                max_batches=1,
            ),
            out_dir=out_dir,
        )
    nu_marcus = result["nu_marcus"]
    mass_marcus = result["mass_marcus"]
    np.testing.assert_allclose((nu_marcus, mass_marcus), (1170.7083, 15.4437), atol=5e5)
