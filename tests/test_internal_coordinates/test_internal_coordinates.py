#!/usr/bin/env python3

import itertools as it
import logging; logging.disable(logging.DEBUG)
import numpy as np
import pytest

from pysisyphus.helpers import geom_from_library
from pysisyphus.InternalCoordinates import RedundantCoords
from pysisyphus.calculators.XTB import XTB
from pysisyphus.optimizers.SteepestDescent import SteepestDescent
from pysisyphus.optimizers.BFGS import BFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.ConjugateGradient import ConjugateGradient

from printer import print_gaussian_ints, read_gaussian_ints, compare_to_gaussian

np.set_printoptions(suppress=True, precision=4)


def get_geom(xyz_fn, coord_type="redund", debug=False):
    geom = geom_from_library(xyz_fn, coord_type=coord_type)
    if debug:
        for pc in geom.internal._prim_internals:
            print(pc.inds+1, pc.val)
    return geom


def get_opt(xyz_fn, coord_type="redund", opt_key="sd"):
    geom = get_geom(xyz_fn, coord_type)
    geom.set_calculator(XTB())
    #return BFGS(geom)
    #return RFOptimizer(geom)
    #return ConjugateGradient(geom)
    opt_dict = {
        "sd": SteepestDescent,
        "cg": ConjugateGradient,
        "rfo": RFOptimizer,
    }
    Opt = opt_dict[opt_key]
    opt_kwargs = {
        "dump": False,
    }
    return Opt(geom, **opt_kwargs)


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
    # opt = get_opt(xyz_fn, opt_key="rfo")
    # opt.dump = True
    # opt.run()

    cart_opt = get_opt(xyz_fn, coord_type="cart", opt_key="rfo")
    cart_opt.dump = True
    cart_opt.run()
    #with open("fe_opt.xyz", "w") as handle:
    #    handle.write(opt.geometry.as_xyz())

    #cart_opt = get_opt(xyz_fn, "cart")
    #cart_opt.run()
    #with open("fe_opt_cart.xyz", "w") as handle:
    #    handle.write(cart_opt.geometry.as_xyz())


def test_biaryl_opt():
    xyz_fn = "01_opt_aus_04_b3lyp_scan_min.xyz"
    opt = get_opt(xyz_fn, opt_key="rfo")
    opt.run()

    cart_opt = get_opt(xyz_fn, coord_type="cart", opt_key="rfo")
    # cart_opt.dump = True
    cart_opt.run()


def test_azetidine():
    """See https://doi.org/10.1063/1.462844 Section. III Examples."""
    xyz_fn = "azetidine_not_opt.xyz"
    geom = get_geom(xyz_fn)
    #import pdb; pdb.set_trace()
    #geom.dihedral_indices
    dihed_inds = geom.internal.dihedral_indices
    dihed_inds = np.concatenate((dihed_inds, ((3,0,8,5), (0,3,5,8), (3,5,8,0))))
    geom.internal.dihedral_indices = dihed_inds
    gaussian_str, pysis_dict = print_gaussian_ints(geom)

    fn = "azetidine_gaussian_internals"
    with open(fn) as handle:
        text = handle.read()
    gaussian_dict = read_gaussian_ints(text)
    #import pdb; pdb.set_trace()
    compare_to_gaussian(pysis_dict, gaussian_dict)
    print("\n".join(gaussian_str))


@pytest.mark.skip(reason="This fails for now ...")
def test_azetidine_opt():
    """See https://doi.org/10.1063/1.462844 Section. III Examples."""
    #xyz_fn = "azetidine_not_opt.xyz"
    #xyz_fn = "azetidine_xtbopt.xyz"
    xyz_fn = "xtbopt_mod2.xyz"
    geom = get_geom(xyz_fn)
    geom = get_geom("/scratch/programme/pysisyphus/xyz_files/azet_bad.xyz")
    dihed_inds = geom.internal.dihedral_indices
    dihed_inds = np.concatenate((dihed_inds, ((3,0,8,5), (0,3,5,8), (3,5,8,0))))
    geom.internal.dihedral_indices = dihed_inds
    geom.internal._prim_internals = geom.internal.calculate(geom._coords)
    print(geom.internal._prim_internals, len(geom.internal._prim_internals))
    print(dihed_inds)
    print(geom.internal.dihedral_indices)
    geom.set_calculator(XTB())
    #return BFGS(geom)
    opt = RFOptimizer(geom)
    #opt = get_opt(xyz_fn)

    #dihed_inds = opt.geometry.internal.dihedral_indices
    #dihed_inds = np.concatenate((dihed_inds, ((3,0,8,5), (0,3,5,8), (3,5,8,0))))
    #opt.geometry.internal.dihedral_indices = dihed_inds

    opt.max_cycles = 100
    opt.dump = True
    opt.run()


def test_runo():
    xyz_fn = "01_run6fbpb_photoprod_guess.xyz"
    geom = get_geom(xyz_fn)
    #import pdb; pdb.set_trace()
    #geom.dihedral_indices
    dihed_inds = geom.internal.dihedral_indices
    dihed_inds = np.concatenate((dihed_inds, ((3,0,8,5), (0,3,5,8), (3,5,8,0))))
    geom.internal.dihedral_indices = dihed_inds
    gaussian_str, pysis_dict = print_gaussian_ints(geom)

    print("\n".join(gaussian_str))


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


def test_co2_linear_opt():
    xyz_fn = "co2_linear.xyz"
    opt = get_opt(xyz_fn)
    opt.run()


def test_ch4():
    xyz_fn = "methane.xyz"
    geom = get_geom(xyz_fn)
    assert_internals(geom, (4, 6, 4))


def test_sf6():
    xyz_fn = "sf6.xyz"
    geom = get_geom(xyz_fn)
    #print(geom.internal.prim_indices)
    gaussian_str, pysis_dict = print_gaussian_ints(geom)
    print("\n".join(gaussian_str))
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


def test_h2o2_opt():
    geom = geom_from_library("h2o2_hf_321g_opt.xyz", coord_type="redund")
    calc = XTB()
    geom.set_calculator(calc)
    H_start = geom.hessian
    ws, vs = np.linalg.eigh(H_start)
    opt_kwargs = {
        "thresh": "gau_tight",
        "hessian_init": "calc",
        # "hessian_recalc": 1,
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    opt.run()
    H_end = geom.hessian
    we, ve = np.linalg.eigh(H_end)
    assert (we > 0).all()


if __name__ == "__main__":
    #test_fluorethylene()
    # test_fluorethylene_opt()
    #test_azetidine()
    #test_azetidine_opt()
    #test_runo()
    #test_h2o()
    #test_h2o_opt()
    #test_h2o_rfopt()
    #test_single_atom_fragments()
    #test_two_fragments()
    #test_hydrogen_bonds()
    #test_co2_linear()
    #test_co2_linear_opt()
    #test_ch4()
    #test_sf6()
    # test_biaryl_opt()
    test_h2o2_opt()
    pass
