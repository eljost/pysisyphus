#!/usr/bin/env python3

import numpy as np

from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.helpers import geom_from_library

np.set_printoptions(suppress=True, precision=4)

# pyscf_ = PySCF(basis="3-21g", pal=4)
pyscf_ = PySCF(basis="3-21g", pal=4)
# pyscf_ = PySCF(basis="3-21g", method="mp2", pal=4)
# geom = geom_from_library("birkholz/vitamin_c.xyz")
# geom = geom_from_library("hcn.xyz")
geom = geom_from_library("hcn_iso_ts.xyz")
geom.set_calculator(pyscf_)
f = geom.forces.reshape(-1, 3)
print("PySCF")
print(f)

H = geom.hessian
print("PySCF hessian")
print(H.reshape(-1, 9))

ref_geom = geom.copy()
from pysisyphus.calculators.Gaussian16 import Gaussian16
g16 = Gaussian16("hf/3-21G", pal=4)
# g16 = Gaussian16("mp2/3-21G", pal=4)
ref_geom.set_calculator(g16)
f_ref = ref_geom.forces.reshape(-1, 3)
print("Gaussian16")
print(f_ref)

# np.testing.assert_allclose(f, f_ref, rtol=5e-3)

H_ref = ref_geom.hessian
print("G16 Hess")
print(H_ref)
