import numpy as np
import pytest

from pysisyphus.config import WF_LIB_DIR
from pysisyphus.wavefunction import Wavefunction


@pytest.mark.parametrize(
    "fn",
    (WF_LIB_DIR / "orca_ch4_def2svp.json",),
)
def test_cgto_eval(fn):
    wf = Wavefunction.from_orca_json(fn)
    n_points = 5
    xyz = np.random.rand(n_points, 3)
    shells = wf.shells
    ref_vals = np.array([shells.eval_single(r) for r in xyz])
    vals = shells.eval(xyz)
    np.testing.assert_allclose(vals, ref_vals.T)
