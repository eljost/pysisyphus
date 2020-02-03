import numpy as np
import pytest

try:
    from pysisyphus.calculators.QCEngine import QCEngine
except ImportError:
    print("QCEngine import failed. Did you install it?")
from pysisyphus.helpers import geom_from_library
from pysisyphus.testing import using_turbomole, using_qcengine, using


@using_turbomole
@using_qcengine
def test_qcengine_turbomole():
    geom = geom_from_library("h2o_guess.xyz")

    qce_kwargs = {
        "program": "turbomole",
        "model": {
            "method": "hf",
            "basis": "def2-SVP",
        },
    }
    qce = QCEngine(**qce_kwargs)

    geom.set_calculator(qce)

    forces = geom.forces
    energy = geom.energy
    norm = np.linalg.norm(forces)

    assert energy == pytest.approx(-75.95615655854)
    assert norm == pytest.approx(0.11354367)


@using("mopac")
@using_qcengine
def test_qcengine_mopac():
    geom = geom_from_library("h2o_guess.xyz")

    qce_kwargs = {
        "program": "mopac",
        "model": {
            "method": "pm6",
            "basis": None,
        },
    }
    qce = QCEngine(**qce_kwargs)

    geom.set_calculator(qce)

    forces = geom.forces
    energy = geom.energy
    print(f"energy={energy:.6f}")
    print(f"forces={forces}")
    norm = np.linalg.norm(forces)

    assert energy == pytest.approx(-0.083446211)
    assert norm == pytest.approx(0.044034367)
