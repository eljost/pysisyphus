import numpy as np
import pytest

from pysisyphus.optimizers.hessian_updates import (
    bfgs_update,
    damped_bfgs_update,
    double_damp,
    sr1_update,
    psb_update,
    flowchart_update,
    mod_flowchart_update,
    bofill_update,
)


@pytest.mark.parametrize(
    "update_func",
    [
        bfgs_update,
        damped_bfgs_update,
        flowchart_update,
        mod_flowchart_update,
        bofill_update,
    ],
)
def test_hessian_updates(update_func):
    N = 3
    dx = np.ones(N)
    dg = np.ones(N)
    H = np.arange(N * N).reshape(-1, N)
    dH = update_func(H, dx, dg)


# double_damp
