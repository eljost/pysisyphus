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
    """Simple test of the Gaussian class. Uses DummyColvar."""
    x0 = 0.0
    g = Gaussian(w=1, s=1, x0=x0)
    val, grad = g.eval(x)
    assert val == pytest.approx(ref_val)
    np.testing.assert_allclose(grad, ref_grad)


def test_gaussian_height0():
    g = Gaussian(w=0.0, s=1, x0=0.0)
    assert g.value(1.0) == pytest.approx(0.0)


"""
Assuming three Gaussians (w=1, s=1), centered at -1, 0 and 1,
the potential at x = 1 will be:
    V = exp(-(1- -1)**2 / 2) + exp(-1**2 / 2) + exp(0)
    V = exp(-2) + exp(-1/2) + 1.0

The gradient is given as
    g = -(1- -1) * exp(-2) + -(1 - 0) * exp(-1/2) + -(1-1) * exp(0)
    g = -2 * exp(-2) - exp(-1/2) + 0
"""
THREE_CENTER_POT = EXP_MINUS2 + EXP_HALF + 1.0
THREE_CENTER_GRAD = -2 * EXP_MINUS2 - EXP_HALF + 0.0


@pytest.mark.parametrize(
    "x, x0, ref_val, ref_grad",
    [
        (
            1.0,
            (-1.0, 0.0, 1.0),
            THREE_CENTER_POT,
            THREE_CENTER_GRAD,
        ),
    ],
)
def test_gaussian_multiple_centers(x, x0, ref_val, ref_grad):
    """Tests Gaussian evaluation with multiple centers."""
    g = Gaussian(w=1, s=1)
    x0 = np.array(x0)
    val = g.value(x, x0)
    assert val == pytest.approx(ref_val)

    grad = g.gradient(x, x0).sum()
    np.testing.assert_allclose(grad, ref_grad)


@pytest.mark.parametrize(
    "x0, ref_val, ref_grad",
    [
        (1, 1.0, np.zeros(6)),
        ((1, 1.1, 1.2), 2.9752111, (0.29554098, 0.0, 0.0, -0.29554098, 0.0, 0.0)),
    ],
)
def test_gaussian_cr_func(x0, ref_val, ref_grad):
    """Test Gaussian evaluation with supplied colvar function."""
    indices = (0, 1)
    c3d = np.zeros((2, 3))
    c3d[0, 0] = 1.0
    cv = CVDistance(indices)
    g = Gaussian(colvar=cv)
    val, grad = g.eval(c3d, x0)

    assert val == pytest.approx(val)
    ref_grad = np.reshape(ref_grad, grad.shape)
    np.testing.assert_allclose(grad, ref_grad)


def test_gaussian_centers():
    g = Gaussian()
    # No Gaussians are present, because x0 == None
    val, grad = g.eval(1)

    assert val == pytest.approx(0)
