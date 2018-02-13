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


if __name__ == "__main__":
    #opt_h2o()
    opt_sf6()
