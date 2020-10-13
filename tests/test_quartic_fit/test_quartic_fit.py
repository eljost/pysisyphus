import numpy as np
import pytest

from pysisyphus.optimizers.poly_fit import quartic_fit


def get_poly_dpoly(deg=4, max_=False):
    factor = -1.0 if max_ else 1.0
    coeffs = factor * np.ones(deg+1)
    poly = np.poly1d(coeffs)
    dpoly = np.polyder(poly)
    return poly, dpoly


def test_quartic_fit_max():
    poly, dpoly = get_poly_dpoly(max_=True)

    dx = 2
    xs = np.linspace(-dx, dx, 100)
    ys = poly(xs)
    true_max = ys.argmax()

    x0 = -1.25
    x1 = 0.25
    e0 = poly(x0)
    e1 = poly(x1)
    g0 = dpoly(x0)
    g1 = dpoly(x1)
    fit_result = quartic_fit(e0, e1, g0, g1, maximize=True)

    # import matplotlib.pyplot as plt
    # fig, ax = plt.subplots(figsize=(8, 8))
    # size = 40
    # ax.plot(xs, ys, label="org.")
    # ax.scatter(xs[true_max], ys[true_max], s=size, c="y", label="true max")
    # ax.scatter((x0, x1), (e0, e1), s=size, c="r", label="calculated")
    # x_fit = (1-fit_result.x) * x0 + fit_result.x * x1
    # ax.scatter(x_fit, fit_result.y, s=size, c="blue", label="maximum")
    # ax.set_xlim(x0-0.5, x1+0.5)
    # ax.set_ylim(0.25*ys.min(), max(1.5*ys.max(), 1))
    # ax.legend()
    # plt.show()

    assert fit_result.x == pytest.approx(0.4973575)
    assert fit_result.y == pytest.approx(-0.79934712)
