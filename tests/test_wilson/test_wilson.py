import numpy as np
import pytest

from pysisyphus.Geometry import Geometry
from pysisyphus.intcoords import (
    Stretch,
    Bend,
    Torsion,
    OutOfPlane,
    LinearBend,
    LinearDisplacement,
    TranslationX,
    TranslationY,
    TranslationZ,
)
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
    q_oop,
    dq_oop,
    d2q_oop,
)
import pysisyphus.intcoords.mp_derivatives as mp_d
from pysisyphus.intcoords.findiffs import fin_diff_prim, fin_diff_B
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
    mp_ref_val = mp_d.q_b(*args)
    mp_ref_grad = mp_d.dq_b(*args)

    assert val == pytest.approx(ref_val)
    assert val == pytest.approx(mp_ref_val)
    np.testing.assert_allclose(grad.flatten(), ref_grad.flatten())
    np.testing.assert_allclose(grad.flatten(), mp_ref_grad.flatten())

    # Code generated 2nd derivative
    dgrad = d2q_b(*args)
    mp_dgrad = mp_d.d2q_b(*args)

    # Finite difference reference values
    ref_dgrad = fin_diff_B(Stretch(indices), coords3d)
    np.testing.assert_allclose(dgrad, ref_dgrad, atol=1e-12)
    np.testing.assert_allclose(mp_dgrad, ref_dgrad, atol=1e-12)


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

    mp_ref_val = mp_d.q_a(*args)
    _mp_ref_grad = mp_d.dq_a(*args)
    mp_ref_grad = np.zeros_like(coords3d)
    mp_ref_grad[indices] = _mp_ref_grad.reshape(-1, 3)

    assert val == pytest.approx(ref_val)
    assert val == pytest.approx(mp_ref_val)
    np.testing.assert_allclose(grad.flatten(), ref_grad.flatten(), atol=1e-12)
    np.testing.assert_allclose(grad.flatten(), mp_ref_grad.flatten(), atol=1e-12)

    # Code generated 2nd derivative
    dgrad = d2q_a(*args)
    mp_dgrad = mp_d.d2q_a(*args)

    # Finite difference reference values
    ref_dgrad = fin_diff_B(Bend(indices), coords3d)
    np.testing.assert_allclose(dgrad, ref_dgrad, atol=1e-9)
    np.testing.assert_allclose(mp_dgrad, ref_dgrad, atol=1e-9)


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

    mp_ref_val = mp_d.q_d(*args)
    _mp_ref_grad = mp_d.dq_d(*args)
    mp_ref_grad = np.zeros_like(coords3d)
    mp_ref_grad[indices] = _mp_ref_grad.reshape(-1, 3)

    # Sign change is not taken into account in q_d
    assert val == pytest.approx(sign * ref_val, abs=1e-8), "Dihedral value"
    assert val == pytest.approx(sign * mp_ref_val, abs=1e-8)
    np.testing.assert_allclose(
        grad.flatten(), sign * ref_grad.flatten(), atol=1e-8, err_msg="1st derivative"
    )
    np.testing.assert_allclose(
        grad.flatten(),
        sign * mp_ref_grad.flatten(),
        atol=1e-8,
        err_msg="1st derivative",
    )

    # Code generated 2nd derivative
    dgrad = sign * d2q_d(*args)
    mp_dgrad = sign * mp_d.d2q_d(*args)

    # Finite difference reference values
    ref_dgrad = fin_diff_B(Torsion(indices), coords3d)
    np.testing.assert_allclose(dgrad, ref_dgrad, atol=1e-8, err_msg="2nd derivative")
    np.testing.assert_allclose(mp_dgrad, ref_dgrad, atol=1e-8, err_msg="2nd derivative")


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
    lb = LinearBend(indices)
    val, grad = lb.calculate(coords3d, gradient=True)
    # Orthogonal direction
    cross_vec = lb.cross_vec
    w = lb._get_orthogonal_direction(coords3d, indices, cross_vec=cross_vec)

    # Reference values, code generated
    args = coords3d[indices].flatten()
    ref_val = q_lb(*args, *w)
    # Reference gradient returned in order [1, 0, 2]
    _ref_grad = dq_lb(*args, *w)
    ref_grad = np.zeros_like(coords3d)
    ref_grad[indices] = _ref_grad.reshape(-1, 3)

    assert val == pytest.approx(ref_val)
    np.testing.assert_allclose(grad.flatten(), ref_grad.flatten(), atol=1e-12)

    # Code generated 2nd derivative
    dgrad = lb.jacobian(coords3d)

    # # Finite difference reference values, only passes for 180°
    # ref_dgrad = fin_diff_B(lb, coords3d)
    # np.testing.assert_allclose(dgrad, ref_dgrad, atol=1e-9)


@pytest.mark.parametrize("dz", np.linspace(-10.0, 10.0, 21))
def test_outofplane(dz):
    indices = [0, 1, 2, 3]

    # Create equilateral triangle
    r = 2
    degs = np.array((120, 240, 360))
    rads = np.deg2rad(degs)
    xs = r * np.cos(rads)
    ys = r * np.sin(rads)
    zs = (0.0, 0.0, 0.0)
    _coords3d = np.stack((xs, ys, zs), axis=1)
    _coords3d -= _coords3d.mean(axis=0)
    coords3d = np.zeros((4, 3))
    coords3d[:3] = _coords3d

    # Add apex atom
    coords3d[3] = (0.0, 0.0, dz)
    atoms = ("C", "C", "C", "C")

    geom = Geometry(atoms, coords3d.flatten())
    # geom.jmol()
    oop = OutOfPlane(indices)
    val, grad = oop.calculate(geom.coords3d, gradient=True)

    # Reference values, code generated
    args = coords3d[indices].flatten()
    ref_val = q_oop(*args)
    assert val == pytest.approx(ref_val)

    # Reference gradient returned in order [0, 1, 2, 3]
    ref_grad = fin_diff_prim(oop, geom.coords3d)
    np.testing.assert_allclose(grad, ref_grad, atol=1e-10)

    # Code generated 2nd derivative
    dgrad = oop.jacobian(coords3d)
    mp_dgrad = mp_d.d2q_oop(*coords3d[indices].flatten())

    # Finite difference reference values
    ref_dgrad = fin_diff_B(oop, coords3d)
    np.testing.assert_allclose(dgrad, ref_dgrad, atol=1e-9)
    np.testing.assert_allclose(mp_dgrad, ref_dgrad, atol=1e-9)


@pytest.mark.parametrize("deg", np.linspace(165, 180, 16))
def test_linear_displacement(deg):
    zmat_str = f"""
    C
    C 1 1.5
    C 1 1.5 2 {deg}
    """
    zmat = zmat_from_str(zmat_str)
    geom = geom_from_zmat(zmat)
    coords3d = geom.coords3d

    indices = [1, 0, 2]
    ld = LinearDisplacement(indices)
    val, grad = ld.calculate(coords3d, gradient=True)

    # First derivative
    ref_row = np.zeros_like(coords3d)
    ref_row[indices] = fin_diff_prim(ld, coords3d).reshape(-1, 3)
    ref_row = ref_row.flatten()
    np.testing.assert_allclose(grad, ref_row, atol=1e-10)

    # Second derivative
    # Code generated 2nd derivative
    dgrad = ld.jacobian(coords3d)
    mp_dgrad = mp_d.d2q_ld(*coords3d[indices].flatten(), *ld.cross_vec)

    # Finite difference reference values
    ref_dgrad = fin_diff_B(ld, coords3d)
    np.testing.assert_allclose(dgrad, ref_dgrad, atol=1e-9)
    np.testing.assert_allclose(mp_dgrad, ref_dgrad, atol=1e-9)


def translation_tester(coords3d, indices):
    for cls in (TranslationX, TranslationY, TranslationZ):
        trans = cls(indices)

        # First derivative
        #   Manually programmed
        val, grad = trans.calculate(coords3d, gradient=True)
        # Finite difference reference values
        ref_grad = np.zeros_like(coords3d)
        ref_grad[indices] = fin_diff_prim(trans, coords3d).reshape(-1, 3)
        np.testing.assert_allclose(grad, ref_grad.flatten())

        # Second derivative
        #   Manually programmed
        dgrad = trans.jacobian(coords3d)
        # Finite difference reference values
        ref_dgrad = fin_diff_B(trans, coords3d)
        np.testing.assert_allclose(dgrad, ref_dgrad)


def test_trans_simple():
    """Simple test for translation coordinates, as proposed by
    Lee-Ping Wang in http://dx.doi.org/10.1063/1.4952956
    """
    zmat_str = f"""
    C
    C 1 3
    """
    zmat = zmat_from_str(zmat_str)
    geom = geom_from_zmat(zmat)
    coords3d = geom.coords3d
    indices = [
        0,
    ]
    translation_tester(coords3d, indices)


def test_trans_complex():
    """See test_trans_simple"""
    hl = 1.5
    dist = 6
    c3d_bot = np.array(((hl, hl, 0.0), (hl, -hl, 0.0), (-hl, -hl, 0.0), (-hl, hl, 0.0)))
    c3d_top = c3d_bot.copy()
    c3d_top[:, 2] += dist
    coords3d = np.concatenate((c3d_bot, c3d_top), axis=0)

    indices_bot = [0, 1, 2, 3]
    translation_tester(coords3d, indices_bot)

    indices_top = [4, 5, 6, 7]
    translation_tester(coords3d, indices_top)
