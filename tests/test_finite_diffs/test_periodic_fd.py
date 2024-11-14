import numpy as np
import pytest

from pysisyphus import finite_diffs as fd


_NUM = 500


@pytest.mark.parametrize(
    "func, ref_func",
    (
        # Test function and its derivative
        (np.cos, lambda x: -np.cos(x)),
        (np.sin, lambda x: -np.sin(x)),
    ),
)
@pytest.mark.parametrize("x_ind", (0, 1, 2, _NUM // 2, _NUM - 4, _NUM - 3, _NUM - 2))
def test_periodic_fd(func, ref_func, x_ind):
    xs, dx = np.linspace(0.0, 2.0 * np.pi, num=_NUM, retstep=True)
    # Drop the last point; this is important to make the grid actually periodic
    xs = xs[:-1]
    ys = func(xs)
    deriv = fd.periodic_fd_2_8(x_ind, ys, dx)
    ref_deriv = ref_func(xs[x_ind])
    assert deriv == pytest.approx(ref_deriv, abs=1e-10)
