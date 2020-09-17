import numpy as np
import pytest

from pysisyphus.intcoords import Stretch, Bend, Torsion, LinearBend
from pysisyphus.intcoords.derivatives import (
    q_b,
    dq_b,
    d2q_b,
    q_a,
    dq_a,
    d2q_a,
    q_d,
    dq_d,
    d2q_d,
    q_lb,
    dq_lb,
    d2q_lb,
)
from pysisyphus.intcoords.findiffs import fin_diff_B
from pysisyphus.io.zmat import geom_from_zmat, zmat_from_str


@pytest.mark.parametrize("length", np.linspace(0.1, 1, 10))
def test_stretch(length):
    indices = [0, 1]

    coords3d = np.array(((0.0, 0.0, 0.0), (0.0, 0.0, length)))
    # Explicitly implemented
    val, grad = Stretch._calculate(coords3d, indices, gradient=True)

    # Reference values, code generated
    args = coords3d[indices].flatten()
    ref_val = q_b(*args)
    ref_grad = dq_b(*args)

    assert val == pytest.approx(ref_val)
    np.testing.assert_allclose(grad.flatten(), ref_grad.flatten())

    # Code generated 2nd derivative
    dgrad = d2q_b(*args)

    # Finite difference reference values
    ref_dgrad = fin_diff_B(Stretch(indices), coords3d)
    np.testing.assert_allclose(dgrad, ref_dgrad, atol=1e-12)


@pytest.mark.parametrize("deg", np.linspace(5, 175, num=35))
def test_bend(deg):
    indices = [1, 0, 2]

    zmat_str = f"""
    C
    C 1 1.5
    C 1 1.5 2 {deg}
    """.strip()
    zmat = zmat_from_str(zmat_str)
    geom = geom_from_zmat(zmat)
    coords3d = geom.coords3d

    # Explicitly implemented
    # Gradient returned in order [0, 1, 2]
    val, grad = Bend._calculate(coords3d, indices, gradient=True)

    # Reference values, code generated
    args = coords3d[indices].flatten()
    ref_val = q_a(*args)
    # Reference gradient returned in order [1, 0, 2]
    _ref_grad = dq_a(*args)
    ref_grad = np.zeros_like(coords3d)
    ref_grad[indices] = _ref_grad.reshape(-1, 3)

    assert val == pytest.approx(ref_val)
    np.testing.assert_allclose(grad.flatten(), ref_grad.flatten(), atol=1e-12)

    # Code generated 2nd derivative
    dgrad = d2q_a(*args)

    # Finite difference reference values
    ref_dgrad = fin_diff_B(Bend(indices), coords3d)
    np.testing.assert_allclose(dgrad, ref_dgrad, atol=1e-9)


@pytest.mark.parametrize(
    "dihedral",
    [
        # First derivative fails alread for 1e-3, -1e-3
        # Fails for ~ 180, ~ 0 and ~ -180
        179,
        140,
        100,
        60,
        20,
        1,
        1,
        -20,
        -60,
        -100,
        -140,
        -179,
    ],
)
def test_torsion(dihedral):
    indices = [3, 2, 0, 1]
    zmat_str = f"""
    C
    C 1 1.
    C 1 1. 2 135.
    C 3 1. 1 135. 2 {dihedral}
    """.strip()
    zmat = zmat_from_str(zmat_str)
    geom = geom_from_zmat(zmat)
    coords3d = geom.coords3d

    # Explicitly implemented
    # Gradient returned in order [3, 2, 0, 1]
    val, grad = Torsion._calculate(coords3d, indices, gradient=True)
    sign = np.sign(val)

    # Reference values, code generated
    args = coords3d[indices].flatten()
    ref_val = q_d(*args)
    # Reference gradient returned in order [3, 2, 0, 1]
    _ref_grad = dq_d(*args)
    ref_grad = np.zeros_like(coords3d)
    ref_grad[indices] = _ref_grad.reshape(-1, 3)

    # Sign change is not taken into account in q_d
    assert val == pytest.approx(sign * ref_val, abs=1e-8), "Dihedral value"
    np.testing.assert_allclose(
        grad.flatten(), sign * ref_grad.flatten(), atol=1e-8, err_msg="1st derivative"
    )

    # Code generated 2nd derivative
    dgrad = sign * d2q_d(*args)

    # Finite difference reference values
    ref_dgrad = fin_diff_B(Torsion(indices), coords3d)
    np.testing.assert_allclose(dgrad, ref_dgrad, atol=1e-8, err_msg="2nd derivative")


@pytest.mark.parametrize("deg", np.linspace(165, 180, num=16))
def test_linear_bend(deg):
    indices = [1, 0, 2]

    zmat_str = f"""
    C
    C 1 1.5
    C 1 1.5 2 {deg}
    """.strip()
    zmat = zmat_from_str(zmat_str)
    geom = geom_from_zmat(zmat)
    coords3d = geom.coords3d

    # Explicitly implemented
    # Gradient returned in order [0, 1, 2]
    val, grad = LinearBend._calculate(coords3d, indices, gradient=True)
    # Orthogonal direction
    w = LinearBend._get_orthogonal_direction(coords3d, indices)

    # Reference values, code generated
    args = coords3d[indices].flatten()
    ref_val = q_lb(*args, *w)
    # Reference gradient returned in order [1, 0, 2]
    _ref_grad = dq_lb(*args, *w)
    ref_grad = np.zeros_like(coords3d)
    ref_grad[indices] = _ref_grad.reshape(-1, 3)

    assert val == pytest.approx(ref_val)
    np.testing.assert_allclose(grad.flatten(), ref_grad.flatten(), atol=1e-12)

    # # Code generated 2nd derivative
    # dgrad = d2q_lb(*args, *w)

    # # Finite difference reference values
    # ref_dgrad = fin_diff_B(LinearBend(indices), coords3d)
    # np.testing.assert_allclose(dgrad, ref_dgrad, atol=1e-9)
