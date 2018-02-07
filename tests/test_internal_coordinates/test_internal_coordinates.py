#!/usr/bin/env python3

import numpy as np

from pysisyphus.helpers import geom_from_library
from pysisyphus.InternalCoordinates import get_B_mat, get_B_inv, backtransform
from pysisyphus.calculators.XTB import XTB


np.set_printoptions(suppress=True, precision=4)


def test_fluorethylene():
    # Fluorethylene, see [2] for geometry
    geom = geom_from_library("fluorethylene.xyz")
    forces_fn = "fe_forces"
    #fe_geom.set_calculator(XTB())
    #forces = fe_geom.forces
    #np.savetxt(forces_fn, forces)
    forces = np.loadtxt(forces_fn)
    print(forces)

    fe_B = get_B_mat(geom)
    B_mat_inv = get_B_inv(fe_B, geom.atoms)
    # Transform forces to internal coordinates
    int_forces = B_mat_inv.dot(forces)
    print(int_forces)
    int_step = 0.5*int_forces
    #max_step = max(abs(step))
    #if max_step > 0.04:
    #        step /= max_step
    backtransform(geom.coords, int_step, B_mat_inv)
    #assert len(fe_inds) == 5
    #assert len(fe_bends) == 6
    #assert len(fe_dihedrals) == 4


def run():
    """
    h2o_geom = geom_from_library("h2o.xyz")
    h2o_B = get_B_mat(h2o_geom)
    #assert len(h2o_inds) == 2
    #assert len(h2o_bends) == 1

    benzene_geom = geom_from_library("benzene_bp86sto3g_opt.xyz")
    benezen_B = get_B_mat(benzene_geom)
    #assert len(benzene_inds) == 12
    """


    """
    # PT H2O
    pt_geom = geom_from_library("h2o_pt.xyz")
    h2o_pt_B = get_B_mat(pt_geom)

    # H2O2, 1 Dihedral
    print("h2o2")
    h2o2_geom = geom_from_library("h2o2_hf_321g_opt.xyz")
    h2o2_B = get_B_mat(h2o2_geom)#, save="h2o2.bmat", tm_format=True)
    """

    """
    # Two fragments
    print("Two fragments")
    two_frags = geom_from_library("h2o2_h2o_fragments.xyz")
    two_frags_B = get_B_mat(two_frags)
    """

if __name__ == "__main__":
    test_fluorethylene()
    #run()
