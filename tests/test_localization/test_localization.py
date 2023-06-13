import h5py
import numpy as np
import pytest

from pysisyphus.config import WF_LIB_DIR
from pysisyphus.linalg import matrix_power
from pysisyphus.wavefunction import (
    cholesky,
    edmiston_ruedenberg,
    foster_boys,
    pipek_mezey,
    Wavefunction,
)


def from_checkpoint(h5_fn):
    with h5py.File(h5_fn, "r") as handle:
        C = handle["Ca"][:]
        S = handle["S"][:]
        int3c2e = handle["3c2e"][:]
        int2c2e = handle["2c2e"][:]
        occs = handle.attrs["occs"]
        assert len(occs) == 1
    Luv = int3c2e.dot(matrix_power(int2c2e, -0.5))
    return C, Luv, S, occs


def test_edmiston_ruedenberg(this_dir):
    h5_fn = this_dir / "ch4_checkpoint.h5"
    C, Luv, S, occs = from_checkpoint(h5_fn)
    nocc = occs[0]
    Luv = Luv.transpose(2, 0, 1)
    Crot = C[:, :nocc]

    # Pre-localize using Cholesky decomposition
    Crot = cholesky(Crot)
    Crot, cost_func = edmiston_ruedenberg(Crot, Luv, S)
    diag = np.diag(Crot.T @ S @ Crot)
    np.testing.assert_allclose(diag, np.ones_like(diag))
    assert cost_func == pytest.approx(6.086131)


def test_pipek_mezey():
    wf = Wavefunction.from_orca_json(WF_LIB_DIR / "orca_n2o4_sto3g.json")
    result = pipek_mezey(wf)
    assert result.is_converged
    assert result.P == pytest.approx(18.59696423)
    assert result.C.shape == (30, 23)


def test_foster_boys():
    wf = Wavefunction.from_orca_json(WF_LIB_DIR / "orca_hcho_def2svp.json")
    result = foster_boys(wf)
    assert result.is_converged
    assert result.P == pytest.approx(277.46881)
    assert result.C.shape == (38, 8)
