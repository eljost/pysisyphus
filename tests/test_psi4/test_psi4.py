#!/usr/bin/env python3

from time import time

from pysisyphus.helpers import geom_from_library
from pysisyphus.calculators import Psi4

geom = geom_from_library("hcn_iso_ts.xyz")

psi4_kwargs = {
    "pal": 4,
    "mem": 2000,
    # "method": "hf",
    "method": "b3lyp",
    "basis": "def2-svp",
    # "mult": 1,
    # "charge": 2,
}
psi4 = Psi4(**psi4_kwargs)
geom.set_calculator(psi4)
print(psi4.base_cmd)
# en = geom.energy
# print(en)
f = geom.forces
print(f)
e = geom.energy
print(e)

start = time()
h = geom.hessian
end = time()
print(h)
dur = end - start
print("hess calc took", int(dur), "seconds")

import numpy as np
det = np.linalg.det(h)
print("det", det)
