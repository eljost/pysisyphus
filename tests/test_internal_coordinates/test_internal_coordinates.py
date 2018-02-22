#!/usr/bin/env python3

import logging; logging.disable(logging.DEBUG)
import numpy as np
import pytest

from pysisyphus.helpers import geom_from_library
from pysisyphus.InternalCoordinates import RedundantCoords, DelocalizedCoords
from pysisyphus.calculators.XTB import XTB
from pysisyphus.optimizers.SteepestDescent import SteepestDescent
from pysisyphus.optimizers.BFGS import BFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.ConjugateGradient import ConjugateGradient


np.set_printoptions(suppress=True, precision=4)


def get_geom(xyz_fn, coord_type="redund", debug=False):
    geom = geom_from_library(xyz_fn, coord_type=coord_type)
    if debug:
        for pc in geom.internal._prim_coords:
            print(pc.inds+1, pc.val)
    return geom


def get_opt(xyz_fn, coord_type="redund"):
    geom = get_geom(xyz_fn, coord_type)
    geom.set_calculator(XTB())
    #return BFGS(geom)
    return RFOptimizer(geom)
    #return ConjugateGradient(geom)
    #return SteepestDescent(geom)


def assert_internals(geom, lengths):
    bonds, bends, dihedrals = lengths
    assert len(geom.internal.bond_indices) == bonds
    assert len(geom.internal.bending_indices) == bends
    assert len(geom.internal.dihedral_indices) == dihedrals


def test_fluorethylene():
    xyz_fn = "fluorethylene.xyz"
    geom = get_geom(xyz_fn)
    assert_internals(geom, (5, 6, 4))


def test_fluorethylene_opt():
    xyz_fn = "fluorethylene.xyz"
    #opt = get_opt(xyz_fn)
    #opt.run()

    cart_opt = get_opt(xyz_fn, coord_type="cart")
    cart_opt.dump = True
    cart_opt.run()
    #with open("fe_opt.xyz", "w") as handle:
    #    handle.write(opt.geometry.as_xyz())

    #cart_opt = get_opt(xyz_fn, "cart")
    #cart_opt.run()
    #with open("fe_opt_cart.xyz", "w") as handle:
    #    handle.write(cart_opt.geometry.as_xyz())


@pytest.mark.skip(reason="This fails for now ...")
def test_azetidine_opt():
    """See https://doi.org/10.1063/1.462844 Section. III Examples."""
    xyz_fn = "azetidine_not_opt.xyz"
    #geom = get_geom(xyz_fn)
    opt = get_opt(xyz_fn)
    opt.dump = True
    geom = opt.geometry
    print("internal coordinates", len(geom.internal._prim_coords))
    opt.run()


def test_h2o():
    xyz_fn = "h2o.xyz"
    geom = get_geom(xyz_fn)
    assert_internals(geom, (2, 1, 0))


def test_h2o_opt():
    xyz_fn = "h2o.xyz"
    opt = get_opt(xyz_fn)
    opt.run()

    cart_opt = get_opt(xyz_fn, "cart")
    cart_opt.run()


def test_h2o_rfopt():
    xyz_fn = "h2o.xyz"
    geom = get_geom(xyz_fn)
    #print(geom.internal.B)
    #print()
    #print(geom.internal.Bt_inv)
    #print()
    #print(geom.internal.P)
    geom.set_calculator(XTB())
    opt = RFOptimizer(geom)
    #import pdb; pdb.set_trace()
    opt.run()


def test_single_atom_fragments():
    xyz_fn = "single_atom_fragments.xyz"
    geom = get_geom(xyz_fn)
    bi = geom.internal.bond_indices.tolist()
    assert ([4, 5] in bi) and ([3, 6] in bi) and ([5, 6] in bi)
    assert(len(geom.internal.fragments) == 3)


def test_two_fragments():
    xyz_fn = "h2o2_h2o_fragments.xyz"
    geom = get_geom(xyz_fn)
    assert(len(geom.internal.fragments) == 2)


def test_hydrogen_bonds():
    xyz_fn = "hydrogen_bond.xyz"
    geom = get_geom(xyz_fn)
    h_bonds = geom.internal.hydrogen_bond_indices.tolist()
    assert ([7, 9] in h_bonds) and ([8, 5] in h_bonds)


def test_co2_linear():
    xyz_fn = "co2_linear.xyz"
    geom = geom_from_library(xyz_fn, coord_type="redund")
    assert_internals(geom, (2, 1, 0))


def test_ch4():
    xyz_fn = "methane.xyz"
    geom = get_geom(xyz_fn)
    assert_internals(geom, (4, 6, 4))


def test_sf6():
    xyz_fn = "sf6.xyz"
    geom = get_geom(xyz_fn)
    assert_internals(geom, (6, 15, 20))


def run():
    """
    benzene_geom = geom_from_library("benzene_bp86sto3g_opt.xyz")
    benezen_B = get_B_mat(benzene_geom)
    #assert len(benzene_inds) == 12

    # PT H2O
    pt_geom = geom_from_library("h2o_pt.xyz")
    h2o_pt_B = get_B_mat(pt_geom)

    # H2O2, 1 Dihedral
    print("h2o2")
    h2o2_geom = geom_from_library("h2o2_hf_321g_opt.xyz")
    h2o2_B = get_B_mat(h2o2_geom)#, save="h2o2.bmat", tm_format=True)
    """


if __name__ == "__main__":
    #test_fluorethylene()
    test_fluorethylene_opt()
    #test_azetidine_opt()
    #test_h2o()
    #test_h2o_opt()
    #test_h2o_rfopt()
    #test_single_atom_fragments()
    #test_two_fragments()
    #test_hydrogen_bonds()
    #test_co2_linear()
    #test_ch4()
    #test_sf6()
    pass
