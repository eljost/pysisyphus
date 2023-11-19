import numpy as np
import pytest

from pysisyphus.linalg import are_collinear


@pytest.mark.parametrize(
    "points, result",
    (
        (((0.0, 0.0, 0.0),), False),
        (((0.0, 0.0, 0.0), (0.0, 0.0, 1.0)), True),
        (((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 0.0, 2.0)), True),
        (((0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 4.0, 0.0)), False),
        (((0.0,), (0.1,)), True),
    ),
)
def test_are_collinear(points, result):
    points = np.array(points)
    assert are_collinear(points) == result
