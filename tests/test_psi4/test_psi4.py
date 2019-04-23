#!/usr/bin/env python3

from pysisyphus.helpers import geom_from_library
from pysisyphus.calculators import Psi4

geom = geom_from_library("hcn_iso_ts.xyz")

psi4_kwargs = {
    "pal": 4,
    "mem": 2000,
    "method": "hf",
    "basis": "sto-3g",
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

h = geom.hessian
print(h)

