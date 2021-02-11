from pysisyphus.dynamics.colvars import CVDistance, CVBend, CVTorsion

import pytest
import numpy as np


@pytest.mark.parametrize(
    "cls, indices, ref_val",
    [
        (CVDistance, (0, 1), 1.0),
        (CVBend, (0, 1, 2), np.pi / 2),
        (CVTorsion, (0, 1, 2, 3), -2.1025204),
    ],
)
def test_prim_internal_cv(cls, indices, ref_val):
    c3d = np.array(
        (
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (1.0, -1.0, 0.0),
            (2.0, -0.5, 1.7),
        )
    )

    # Check if autograd gradient is used
    cv = cls(indices)
    agrad_cv = cls(indices, force_agrad=True)
    assert agrad_cv.agrad

    val, grad = cv.eval(c3d)
    agrad_grad = agrad_cv.gradient(c3d)

    assert val == pytest.approx(ref_val)
    np.testing.assert_allclose(grad, agrad_grad, atol=1e-12)


def test_cvbend():
    c3d = np.array(
        (
            (-1.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
        )
    )
    indices = (0, 1, 2)
    cvbend = CVBend(indices)
    with pytest.raises(ZeroDivisionError):
        val, grad = cvbend.eval(c3d)
