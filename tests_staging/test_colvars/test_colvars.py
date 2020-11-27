from pysisyphus.dynamics.colvars import CVDistance

import pytest
import numpy as np
import autograd
import autograd.numpy as anp


def test_cvdistance():
    coords = np.array((
        (1.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
    ))
    indices = (0, 1)
    cv = CVDistance(indices)
    assert cv.agrad == False
    agrad_cv = CVDistance(indices, force_agrad=True)
    assert agrad_cv.agrad

    ref_val = 1.0
    val, grad = cv.eval(coords)
    agrad_val, agrad_grad = agrad_cv.eval(coords)

    assert val == pytest.approx(ref_val)
    np.testing.assert_allclose(grad, agrad_grad)#coords)
