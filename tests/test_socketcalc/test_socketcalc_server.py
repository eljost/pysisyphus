import numpy as np
import pytest

from pysisyphus.helpers import geom_loader
from pysisyphus.calculators.SocketCalc import SocketCalc

@pytest.mark.skip
def test_socketcalc():
    geom = geom_loader("h2o.xyz")
    calc = SocketCalc()
    geom.set_calculator(calc)
    energy = geom.energy

    # Assuming XTB-GFN2
    assert energy == pytest.approx(-5.070431327395)
    print("Got energy", energy)

    forces = geom.forces
    norm_forces = np.linalg.norm(forces)
    assert norm_forces == pytest.approx(0.00578811)
    print("Got forces", forces)

    hessian = geom.hessian
    print("Got hessian", hessian)


if __name__ == "__main__":
    test_socketcalc()
