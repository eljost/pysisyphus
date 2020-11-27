from pysisyphus.dynamics.colvars import CVDistance, CVBend

import pytest
import numpy as np
import autograd
import autograd.numpy as anp


# def test_cvdistance():
@pytest.mark.parametrize(
    "cls, indices, ref_val", [
        (CVDistance, (0, 1), 1.0),
        (CVBend, (0, 1, 2), np.pi / 2),
    ]
)
def test_prim_internal_cv(cls, indices, ref_val):
    coords = np.array((
        (0.0,  0.0, 0.0),
        (1.0,  0.0, 0.0),
        (1.0, -1.0, 0.0),
        (2.0, -1.0, 0.0),
    ))

    # Check if autograd gradient is used
    cv = cls(indices)
    assert cv.agrad == False
    agrad_cv = cls(indices, force_agrad=True)
    assert agrad_cv.agrad

    val, grad = cv.eval(coords)
    agrad_grad = agrad_cv.gradient(coords)

    assert val == pytest.approx(ref_val)
    np.testing.assert_allclose(grad, agrad_grad)
