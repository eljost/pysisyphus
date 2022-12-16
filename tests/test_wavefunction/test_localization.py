import pytest

from pysisyphus.config import WF_LIB_DIR
from pysisyphus.wavefunction import foster_boys, pipek_mezey, Wavefunction


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
