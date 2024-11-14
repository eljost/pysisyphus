import numpy as np
import pytest

from pysisyphus.constants import BOHR2ANG
from pysisyphus.helpers import geom_loader
from pysisyphus.io import geom_from_hessian

from pysisyphus.franckcondon import (
    unitless_displs_from_eigensystem,
    duschinsky,
    DuschinskyRef,
    get_axis_switch,
)


@pytest.mark.parametrize(
    "ref, K0_ref, K_unitless_ref",
    (
        (DuschinskyRef.INITIAL, 0.34956251, 1.176230738),
        (DuschinskyRef.FINAL, -0.34956251, -1.176230738),
    ),
)
def test_duschinsky(ref, K0_ref, K_unitless_ref):
    """Test data taken from strawberry fields program."""
    L_init = np.array([[-0.28933191], [0.0], [0.0], [0.95711104], [0.0], [0.0]])
    L_final = np.array([[-0.28933191], [0.0], [0.0], [0.95711104], [0.0], [0.0]])
    coords3d_init = (
        np.array([-0.0236, 0.0, 0.0, 1.2236, 0.0, 0.0]).reshape(-1, 3) / BOHR2ANG
    )
    coords3d_final = (
        np.array([0.0, 0.0, 0.0, 1.4397, 0.0, 0.0]).reshape(-1, 3) / BOHR2ANG
    )
    masses = np.array((11.0093, 1.0078))
    res = duschinsky(
        L_init,
        coords3d_init,
        L_final,
        coords3d_final,
        masses=masses,
        reference=ref,
        with_axis_switch=False,
    )
    K = res.K
    assert K[0] == pytest.approx(K0_ref)
    wavenums = np.array((1363.2,))
    K_unitless = res.to_K_unitless(wavenums)
    assert K_unitless[0] == pytest.approx(K_unitless_ref)


def test_axis_switching():
    xyz_neutral = """
    9

    Co 0.000004 -0.000318 -0.212456
    C -0.000015 -0.000402 1.628417
    C -1.555843 -0.897775 -0.516878
    C 0.000098 1.795842 -0.516414
    C 1.555765 -0.897933 -0.516847
    O -0.000042 -0.000817 2.766589
    O -2.530524 -1.459219 -0.703171
    O 0.000157 2.920787 -0.701967
    O 2.530393 -1.459475 -0.70312
    """
    xyz_anion = """
    9

    Co 0 0 0
    C 1.01298 1.01298 1.01298
    C -1.01298 -1.01298 1.01298
    C -1.01298 1.01298 -1.01298
    C 1.01298 -1.01298 -1.01298
    O 1.683685 1.683685 1.683685
    O -1.683685 -1.683685 1.683685
    O -1.683685 1.683685 -1.683685
    O 1.683685 -1.683685 -1.683685
    """
    geom_neutral = geom_loader(xyz_neutral)
    geom_anion = geom_loader(xyz_anion)

    # From neutral to anion/ Neutral is initial state
    coords3d_init = geom_neutral.coords3d
    coords3d_final = geom_anion.coords3d
    masses = geom_neutral.masses

    axis_switch = get_axis_switch(coords3d_init, coords3d_final, masses)
    T0 = axis_switch.T0
    # Table 6 from https://pubs.acs.org/doi/pdf/10.1021/jp004230b
    T0_ref = np.array(
        (
            (0.70708, -0.40846, 0.57723),
            (0.00003, 0.81632, 0.5776),
            (-0.70713, -0.40839, 0.57722),
        ),
    )
    np.testing.assert_allclose(T0.T, T0_ref, atol=5e-6)

    # Determine again ... should give the unit matrix now
    axis_switch = get_axis_switch(axis_switch.coords3d_init_rot, coords3d_final, masses)
    T0 = axis_switch.T0
    np.testing.assert_allclose(T0.T, np.eye(3), atol=1e-14)


def test_co_duschinsky(this_dir):
    geom_neutral = geom_from_hessian(this_dir / "coco4_neutral.h5")
    geom_anion = geom_from_hessian(this_dir / "coco4_anion.h5")

    # From neutral to anion/ Neutral is initial state
    geom_init = geom_neutral
    geom_final = geom_anion

    res = duschinsky(geom_init, geom_final)

    print()
    print(f"det(J)={np.linalg.det(res.J):.8f}")
    nus_init, _, L_init, _ = geom_init.get_normal_modes()
    for i, nu in enumerate(nus_init):
        print(f"ν_{i+1:02d} {nu: >12.4f} {res.K[i]: >10.6f}")
    assert abs(res.K[4]) == pytest.approx(2.773520)
    assert abs(res.K[9]) == pytest.approx(5.916038)


def test_unitless_displs_from_eigensystem():
    n = 2
    mw_gradient = np.arange(n)
    H = np.eye(n)
    w = np.diag(H)
    v = H.copy()
    displs = unitless_displs_from_eigensystem(mw_gradient, w, v)
    displs_ref = np.array((0.0, -6.53416392))
    np.testing.assert_allclose(displs, displs_ref)
