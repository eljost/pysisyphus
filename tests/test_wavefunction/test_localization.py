from pysisyphus.config import WF_LIB_DIR
from pysisyphus.wavefunction import pipek_mezey, Wavefunction


def test_pipek_mezey():
    wf = Wavefunction.from_orca_json(WF_LIB_DIR / "orca_n2o4_sto3g.json")
    Ca_loc = pipek_mezey(wf)
    assert Ca_loc.shape == (30, 23)
