#!/usr/bin/env python3

import logging; logging.disable(logging.DEBUG)
import numpy as np

from pysisyphus.helpers import geom_from_library
from pysisyphus.InternalCoordinates import RedundantCoords, DelocalizedCoords
from pysisyphus.calculators.XTB import XTB
from pysisyphus.optimizers.SteepestDescent import SteepestDescent
from pysisyphus.optimizers.BFGS import BFGS
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


np.set_printoptions(suppress=True, precision=4)


def run_opt(xyz_fn, coord_type="redund"):
    geom = geom_from_library(xyz_fn, coord_type=coord_type)
    geom.set_calculator(XTB())
    opt = BFGS(geom)
    #opt = RFOptimizer(geom)
    opt.run()

    return geom


def test_fluorethylene():
    xyz_fn = "fluorethylene.xyz"
    geom = run_opt(xyz_fn)

    geom = run_opt(xyz_fn, "cart")

    #assert len(fe_inds) == 5
    #assert len(fe_bends) == 6
    #assert len(fe_dihedrals) == 4


def test_h2o():
    xyz_fn = "h2o.xyz"
    geom = run_opt(xyz_fn)

    geom = run_opt(xyz_fn, "cart")

    #assert len(h2o_inds) == 2
    #assert len(h2o_bends) == 1


def test_single_atom_fragments():
    xyz_fn = "single_atom_fragments.xyz"
    geom = geom_from_library(xyz_fn, coord_type="redund")
    bi = geom.internal.bond_indices
    assert ([4, 5] in bi) and ([3, 6] in bi) and ([5, 6] in bi)


def test_sf6():
    xyz_fn = "sf6.xyz"
    geom = geom_from_library(xyz_fn, coord_type="redund")
    #import pdb; pdb.set_trace()
    for pc in geom.internal._prim_coords:
        print(pc.inds, pc.val)


def test_hydrogen_bonds():
    xyz_fn = "hydrogen_bond.xyz"
    geom = geom_from_library(xyz_fn, coord_type="redund")
    h_bonds = geom.internal.hydrogen_bond_indices
    assert(len(h_bonds) == 2)


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

    # Two fragments
    print("Two fragments")
    two_frags = geom_from_library("h2o2_h2o_fragments.xyz")
    """


if __name__ == "__main__":
    #test_fluorethylene()
    #test_h2o()
    #test_single_atom_fragments()
    #test_sf6()
    test_hydrogen_bonds()
    #run()
