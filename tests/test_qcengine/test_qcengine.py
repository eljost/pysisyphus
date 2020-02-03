import numpy as np
import pytest

from pysisyphus.calculators.QCEngine import QCEngine
from pysisyphus.helpers import geom_from_library
from pysisyphus.testing import using_turbomole, using_qcengine


@using_turbomole
@using_qcengine
def test_qcengine():
    geom = geom_from_library("h2o_guess.xyz")

    qce_kwargs = {
        "program": "turbomole",
        "model": {
            "method": "hf",
            "basis": "def2-SVP",
        },
        "keywords": {},
    }
    qce = QCEngine(**qce_kwargs)

    geom.set_calculator(qce)

    forces = geom.forces
    energy = geom.energy
    print(f"energy={energy:.6f}")
    print(f"forces={forces}")
    norm = np.linalg.norm(forces)

    assert energy == pytest.approx(-75.95615655854)
    assert norm == pytest.approx(0.11354367)
