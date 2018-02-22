#!/usr/bin/env python3

import numpy as np

np.set_printoptions(suppress=True, precision=4)
from pyscf.geomopt.berny_solver import optimize

from pyscf import gto, scf


def opt_h2o():
    mol = gto.M(atom="""
     O                 -0.00000000   -0.11081188    0.00000000
     H                  0.78397589    0.44324751    0.00000000
     H                 -0.78397589    0.44324751    0.00000000
     """, basis="sto-3g")

    mf = scf.RHF(mol)
    mol1 = optimize(mf, maxsteps=1)


def opt_sf6():
    sf6_mol = gto.M(atom="""
     S                 -0.59356137    0.59356136    0.00000000
     F                 -0.59356137    0.59356136   -1.59000000
     F                 -0.59356137    2.18356136    0.00000000
     F                 -2.18356137    0.59356136    0.00000000
     F                 -0.59356137    0.59356136    1.59000000
     F                 -0.59356137   -0.99643864    0.00000000
     F                  0.99643863    0.59356136    0.00000000
    """, basis="sto-3g")
    mf = scf.RHF(sf6_mol)
    mol1 = optimize(mf)#, maxsteps=1)


def opt_azetidine():
    mol = gto.M(atom="""
     C                 -2.34863259    0.70994422    0.14975543
     H                 -1.82664215   -0.18774343   -0.10828112
     H                 -3.39537825    0.49299009    0.10333942
     N                 -1.96488182    1.28417466    1.38787045
     H                 -0.98378042    1.23964802    1.57617181
     C                 -2.38604561    2.48589673    0.76463760
     H                 -3.44072555    2.64557427    0.84862153
     H                 -1.90158885    3.36988967    1.12346645
     C                 -2.01779803    2.00680465   -0.70250775
     H                 -2.49979795    1.48773285   -1.50446820
     H                 -2.52937959    2.89193912   -1.01829449
     """, basis="sto-3g")
    mf = scf.RHF(mol)
    mol1 = optimize(mf)


def opt_fluorethylene():
    mol = gto.M(atom="""
        C       -0.061684   0.673790    0.000000
        C       -0.061684  -0.726210    0.000000
        F        1.174443   1.331050    0.000000
        H       -0.927709   1.173790    0.000000
        H       -0.927709  -1.226210    0.000000
        H        0.804342  -1.226210    0.000000
    """, basis="sto-3g")
    mf = scf.RHF(mol)
    mol1 = optimize(mf)


def opt_no():
    mol = gto.M(atom="""
        N       0.0     0.0     0.0
        O       0.0     0.0     1.5
    """, basis="sto-3g", spin=1)
    mf = scf.UKS(mol)
    mol1 = optimize(mf)

if __name__ == "__main__":
    #opt_h2o()
    #opt_sf6()
    #opt_azetidine()
    #opt_fluorethylene()
    opt_no()
