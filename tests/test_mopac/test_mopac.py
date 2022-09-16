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
    # assert energy == pytest.approx(-0.091702313)
    assert energy == pytest.approx(-0.09170261366094087)
    assert norm == pytest.approx(0.0407509836626757)

    nus, *_ = geom.get_normal_modes()
    ref = [1370.787877, 2780.770355, 2786.987219]
    np.testing.assert_allclose(nus, ref, atol=1e-2)
