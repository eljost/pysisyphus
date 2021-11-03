import numpy as np
import pytest

from pysisyphus.calculators import MOPAC
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using


@using("mopac")
def test_mopac():
    geom = geom_loader("lib:h2o.xyz")
    calc = MOPAC()
    geom.set_calculator(calc)
    forces = geom.forces
    norm = np.linalg.norm(forces)
    energy = geom.energy
    assert energy == pytest.approx(-0.091702313)
    np.testing.assert_allclose(norm, 0.0407509836626757)

    nus, *_ = geom.get_normal_modes()
    ref = (1369.84, 2778.85, 2785.06)
    np.testing.assert_allclose(nus, ref, atol=1e-2)
