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


"""
Assuming three Gaussians (w=1, s=1), centered at -1, 0 and 1,
the potential at x = 1 will be:
    V = exp(-(1- -1)**2 / 2) + exp(-1**2 / 2) + exp(0)
    V = exp(-2) + exp(-1/2) + 1.0

The gradient is given as
    g = -(1- -1) * exp(-2) + -(1 - 0) * exp(-1/2) + -(1-1) * exp(0)
    g = -2 * exp(-2) - exp(-1/2) + 0

If we extend this to 2D, with same centers and same x the potential
is just doubled, and the gradient is the same in both dimensions.
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
        (
            (1.0, 1.0),
            (-1.0, 0.0, 1.0, -1.0, 0.0, 1.0),
            2 * THREE_CENTER_POT,
            (THREE_CENTER_GRAD, THREE_CENTER_GRAD),
        ),
    ],
)
def test_gaussian_multiple_centers(x, x0, ref_val, ref_grad):
    g = Gaussian(w=1, s=1)
    x = np.array(x)
    x0 = np.array(x0)
    val = g.value(x, x0)
    assert val == pytest.approx(ref_val)

    grad = g.gradient(x, x0)
    np.testing.assert_allclose(grad, ref_grad)


# @pytest.mark.xfail
def test_gaussian_cr_func():
    indices = (0, 1)
    c3d = np.zeros((2, 3))
    c3d[0, 0] = 1.0
    cv = CVDistance(indices)
    x0 = (1, 1.1, 1.2)
    x = cv.value(c3d)
    g = Gaussian(colvar=cv)
    val, grad = g.eval(c3d, x0)
