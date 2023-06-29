import numpy as np
import pytest


from pysisyphus.drivers.boltzmann import boltzmann_weights


@pytest.mark.parametrize(
    "energies, ref_weights",
    (
        (np.arange(3), (1.0, 0.0, 0.0)),
        (np.ones(4), (0.25, 0.25, 0.25, 0.25)),
        (np.arange(3) * 1e-3, (0.68166, 0.236374, 0.081966)),
    ),
)
def test_boltzmann_weights(energies, ref_weights):
    T = 298.15
    weights = boltzmann_weights(energies, T)
    np.testing.assert_allclose(weights, ref_weights, atol=3e-7)

    assert weights.sum() == pytest.approx(1.0)
