from pysisyphus.config import WF_LIB_DIR
from pysisyphus.wavefunction import foster_boys, pipek_mezey, Wavefunction


def test_pipek_mezey():
    wf = Wavefunction.from_orca_json(WF_LIB_DIR / "orca_n2o4_sto3g.json")
    Ca_loc = pipek_mezey(wf)
    assert Ca_loc.shape == (30, 23)


def test_foster_boys():
    wf = Wavefunction.from_orca_json(WF_LIB_DIR / "orca_hcho_def2svp.json")
    Ca_loc = foster_boys(wf)
    assert Ca_loc.shape == (38, 8)
