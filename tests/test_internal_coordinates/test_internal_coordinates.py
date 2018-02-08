#!/usr/bin/env python3

import numpy as np

from pysisyphus.helpers import geom_from_library
from pysisyphus.InternalCoordinates import RedundantCoords, DelocalizedCoords
from pysisyphus.calculators.XTB import XTB
from pysisyphus.optimizers.SteepestDescent import SteepestDescent


np.set_printoptions(suppress=True, precision=4)


def base(xyz_fn):
    geom = geom_from_library(xyz_fn)
    geom.set_calculator(XTB())
    rc = RedundantCoords(geom)

    return geom, rc


def test_fluorethylene():
    # Fluorethylene, see [2] for geometry
    xyz_fn = "fluorethylene.xyz"
    geom, rc = base(xyz_fn)
    forces_fn = "fe_forces"
    forces = np.loadtxt(forces_fn)
    #print(forces)
    step = rc.B_inv.dot(forces)
    rc.transform_int_step(step)
    init_hess = rc.get_initial_hessian()
    #print("rho")
    #print(rc.rho)
    #print("guess_hessian")
    #print(init_hess)

    #assert len(fe_inds) == 5
    #assert len(fe_bends) == 6
    #assert len(fe_dihedrals) == 4


def test_h2o():
    xyz_fn = "h2o.xyz"
    geom, rc = base(xyz_fn)
    forces_fn = "h2o_forces"
    #forces = geom.forces
    forces = np.loadtxt(forces_fn)
    #print(forces)
    step = rc.B_inv.dot(forces)
    rc.transform_int_step(step)
    init_hess = rc.get_initial_hessian()
    #print("rho")
    #print(rc.rho)
    #print("guess_hessian")
    #print(init_hess)

    #assert len(h2o_inds) == 2
    #assert len(h2o_bends) == 1


def test_internal_h2o():
    xyz_fn = "h2o.xyz"
    geom = geom_from_library(xyz_fn)
    rc = RedundantCoords(geom.atoms, geom.coords)
    rc.set_calculator(XTB())
    #forces_fn = "h2o_forces"
    #forces = np.loadtxt(forces_fn)
    opt = SteepestDescent(rc)
    opt.run()
    #step = rc.B_inv.dot(forces)
    #rc.transform_int_step(step)
    #init_hess = rc.get_initial_hessian()


def test_h2o_opt():
    geom = geom_from_library("h2o.xyz")
    geom.set_calculator(XTB())
    ic = RedundantCoords(geom)
    for i in range(25):
        forces = geom.forces
        forces_norm = np.linalg.norm(forces)
        print(f"norm(forces) = {forces_norm:1.4f}")
        if forces_norm < 1e-3:
            break
        internal_step = ic.B_inv.dot(forces)
        ic.transform(internal_step)


def run():
    """

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

def test_two_fragments():
    # Two fragments
    print("Two fragments")
    two_frags = geom_from_library("h2o2_h2o_fragments.xyz")
    two_frags_B = get_B_mat(two_frags)

if __name__ == "__main__":
    #test_fluorethylene()
    #test_h2o()
    test_internal_h2o()
    #test_h2o_opt()
    #test_two_fragments()
    #run()
