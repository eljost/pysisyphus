import matplotlib.pyplot as plt
import numpy as np
import pytest

from pysisyphus.config import WF_LIB_DIR
from pysisyphus.constants import AU2MDYNEPERANG
from pysisyphus.drivers.local_force_constants import (
    get_force_constants_from_complice_mat,
    get_local_force_constants,
    local_mode_overlaps,
    stretch_constants,
)
from pysisyphus.Geometry import Geometry
from pysisyphus.intcoords.PrimTypes import Bonds, Bends, PrimTypes as PT
from pysisyphus.io.fchk import geom_with_hessian_from_fchk
from pysisyphus.io.hessian import geom_from_hessian


@pytest.fixture
def water_geom():
    """Psi4 CCSD(T)/aug-cc-pvtz, opt & Hessian"""
    atoms = ("O", "H", "H")
    coords3d = np.array(
        [
            [0.0, 0.0, -0.12359836],
            [0.0, -1.42992089, 0.9807978],
            [-0.0, 1.42992089, 0.9807978],
        ]
    )
    hessian = np.array(
        [
            [
                -7.78469112e-05,
                0.00000000e00,
                0.00000000e00,
                3.89234556e-05,
                0.00000000e00,
                0.00000000e00,
                3.89234556e-05,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                7.02217124e-01,
                0.00000000e00,
                0.00000000e00,
                -3.51108562e-01,
                2.71211283e-01,
                0.00000000e00,
                -3.51108562e-01,
                -2.71211283e-01,
            ],
            [
                0.00000000e00,
                0.00000000e00,
                4.68159823e-01,
                0.00000000e00,
                2.09620950e-01,
                -2.34079911e-01,
                0.00000000e00,
                -2.09620950e-01,
                -2.34079911e-01,
            ],
            [
                3.89234556e-05,
                0.00000000e00,
                0.00000000e00,
                2.14490726e-06,
                0.00000000e00,
                0.00000000e00,
                -4.10683629e-05,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                -3.51108562e-01,
                2.09620950e-01,
                0.00000000e00,
                3.83142772e-01,
                -2.40416117e-01,
                0.00000000e00,
                -3.20342096e-02,
                3.07951667e-02,
            ],
            [
                0.00000000e00,
                2.71211283e-01,
                -2.34079911e-01,
                0.00000000e00,
                -2.40416117e-01,
                2.21793525e-01,
                0.00000000e00,
                -3.07951667e-02,
                1.22863869e-02,
            ],
            [
                3.89234556e-05,
                0.00000000e00,
                0.00000000e00,
                -4.10683629e-05,
                0.00000000e00,
                0.00000000e00,
                2.14490726e-06,
                0.00000000e00,
                0.00000000e00,
            ],
            [
                0.00000000e00,
                -3.51108562e-01,
                -2.09620950e-01,
                0.00000000e00,
                -3.20342096e-02,
                -3.07951667e-02,
                0.00000000e00,
                3.83142772e-01,
                2.40416117e-01,
            ],
            [
                0.00000000e00,
                -2.71211283e-01,
                -2.34079911e-01,
                0.00000000e00,
                3.07951667e-02,
                1.22863869e-02,
                0.00000000e00,
                2.40416117e-01,
                2.21793525e-01,
            ],
        ]
    )
    geom = Geometry(atoms, coords3d)
    geom.cart_hessian = hessian
    return geom


def test_local_force_constants(water_geom):
    geom = water_geom
    geom_redund = geom.copy_all(coord_type="redund")
    B = geom_redund.internal.B
    nus, *_, L = geom.get_normal_modes()
    force_constants_compl, _ = get_force_constants_from_complice_mat(geom.hessian, B)
    force_constants_local, local_modes = get_local_force_constants(geom.hessian, B, L)
    np.testing.assert_allclose(force_constants_compl, force_constants_local, atol=1e-4)
    S, C = local_mode_overlaps(geom.hessian, L, local_modes)
    # fig, ax = plt.subplots()
    # nuss = [f"{nu:4.2f}" for nu in nus]
    # for i, local_mode in enumerate(C):
    # ax.bar(nuss, local_mode, bottom=C[:i].sum(axis=0))
    # plt.show()
    np.testing.assert_allclose(C.sum(axis=0), np.ones(len(nus)))

    typed_prims = geom_redund.internal.typed_prims
    print()
    # TODO: Below, some code corrects the force constant for bends, but it is hardcoded
    # for systems with 3 internal coordinates, 2 bonds and 1 bend.
    valid_prim_types = Bonds + Bends
    for i, ki in enumerate(force_constants_local):
        prim_type, *prim_inds = typed_prims[i]
        if not (prim_type in valid_prim_types):
            continue

        kcomp = force_constants_compl[i]
        if prim_type in Bends:
            bond_1_inds = prim_inds[:2]
            bond_2_inds = prim_inds[1:]
            ind1 = geom_redund.internal.get_index_of_typed_prim((PT.BOND, *bond_1_inds))
            ind2 = geom_redund.internal.get_index_of_typed_prim((PT.BOND, *bond_2_inds))
            bond_length_1 = geom_redund.coords[ind1]
            bond_length_2 = geom_redund.coords[ind2]
            prod = np.product(1 / (bond_length_1 * bond_length_2))
            ki *= prod
            kcomp *= prod
        print(
            f"{i:03d}: ({prim_type}, {prim_inds}), ki={ki*AU2MDYNEPERANG}, "
            f"compliance k={kcomp*AU2MDYNEPERANG}"
        )


def test_bf3(this_dir):
    geom = geom_from_hessian(this_dir / "bf3_ccsd_t_ccpvdz.h5")

    tps = (
        (PT.BOND, 0, 1),
        (PT.BOND, 0, 2),
        (PT.BOND, 0, 3),
        (PT.BEND, 1, 0, 2),
        (PT.BEND, 2, 0, 3),
        (PT.IMPROPER_DIHEDRAL, 1, 0, 2, 3),
    )
    geom_redund = geom.copy_all(
        coord_type="redund",
        coord_kwargs={
            "typed_prims": tps,
        },
    )
    B = geom_redund.internal.B
    nus, *_, L = geom.get_normal_modes()
    force_constants_compl, _ = get_force_constants_from_complice_mat(geom.hessian, B)
    force_constants_local, local_modes = get_local_force_constants(geom.hessian, B, L)
    # np.testing.assert_allclose(force_constants_compl, force_constants_local, atol=1e-4)
    S, C = local_mode_overlaps(geom.hessian, L, local_modes)
    # fig, ax = plt.subplots()
    # nuss = [f"{nu:4.2f}" for nu in nus]
    # for i, local_mode in enumerate(C):
    # ax.bar(nuss, local_mode, bottom=C[:i].sum(axis=0), label=tps[i])
    # ax.legend()
    # plt.show()


@pytest.mark.parametrize(
    "fchk_fn, ref_consts",
    (
        (
            "01_ch4_freq.fchk",
            {
                (PT.BOND, 0, 1): 0.3445713173236299,
                (PT.BOND, 0, 2): 0.3445724116263941,
                (PT.BOND, 0, 3): 0.34457241135707334,
                (PT.BOND, 0, 4): 0.3445693463705333,
            },
        ),
        (
            "03_ch3f_freq.fchk",
            {
                (PT.BOND, 0, 1): 0.34626240686928556,
                (PT.BOND, 0, 2): 0.33262918333434555,
                (PT.BOND, 0, 3): 0.33262952711888816,
                (PT.BOND, 0, 4): 0.3326307722154276,
            },
        ),
        (
            "05_acn_freq.fchk",
            {
                (PT.BOND, 0, 2): 1.2735519829437154,
                (PT.BOND, 1, 2): 0.32808724468158307,
                (PT.BOND, 1, 3): 0.34017220767653955,
                (PT.BOND, 1, 4): 0.34019072729400324,
                (PT.BOND, 1, 5): 0.3401894660558234,
            },
        ),
        (
            "06_acn_non_statpt.fchk",
            {
                (PT.BOND, 0, 2): 1.2743580248319792,
                (PT.BOND, 1, 2): 0.33693854795027894,
                (PT.BOND, 1, 3): 0.3408064569080432,
                (PT.BOND, 1, 4): 0.34061848807477624,
                (PT.BOND, 1, 5): 0.34079552137470376,
            },
        ),
    ),
)
def test_stretch_constants(fchk_fn, ref_consts):
    """Tests reproduce values report in 10.1063/1.4747339 ."""
    geom = geom_with_hessian_from_fchk(WF_LIB_DIR / fchk_fn)
    results = stretch_constants(geom)
    # print(repr(results))
    for rk, rv in ref_consts.items():
        assert results[rk] == pytest.approx(rv)
