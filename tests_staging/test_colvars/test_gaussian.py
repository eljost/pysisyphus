import math

import numpy as np
import pytest

from pysisyphus.dynamics.Gaussian import Gaussian
from pysisyphus.dynamics.colvars import CVDistance


EXP_HALF = math.exp(-1 / 2)
EXP_MINUS2 = math.exp(-2)


@pytest.mark.parametrize(
    "x, ref_val, ref_grad",
    [
        (-1, EXP_HALF, EXP_HALF),
        (0, 1, 0.0),
        (1, EXP_HALF, -EXP_HALF),
    ],
)
def test_gaussian_value(x, ref_val, ref_grad):
    x0 = 0.0
    g = Gaussian(w=1, s=1, x0=x0)
    val, grad = g.eval(x)
    assert val == pytest.approx(ref_val)
    np.testing.assert_allclose(grad, ref_grad)


def test_gaussian_height0():
    g = Gaussian(w=0.0, s=1, x0=0.0)
    assert g.value(1.0) == pytest.approx(0.0)


def test_gaussian_multiple_centers():
    g = Gaussian(w=1, s=1)
    x = 1.0
    x0 = np.array(
        (
            (-1.0,),  # 2 units to the left
            (0.0,),  # 1 unit to the left
            (1.0,),  # directly at x0
        )
    )
    val = g.value(x, x0)
    assert val == pytest.approx(EXP_MINUS2 + EXP_HALF + 1.0)

    grad = g.gradient(x, x0)
    np.testing.assert_allclose(grad, -2 * EXP_MINUS2 + -EXP_HALF + 0.0)
