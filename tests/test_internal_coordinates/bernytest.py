#!/usr/bin/env

import numpy as np

from berny.geomlib import load
from berny.coords import InternalCoords


np.set_printoptions(suppress=True, precision=4)

fn = "/scratch/programme/pysisyphus/xyz_files/h2o.xyz"
with open(fn) as handle:
    mol = load(handle, fmt="xyz")
ic = InternalCoords(mol)
print(mol)
print(ic)
bmat = ic.B_matrix(mol)
print(bmat)

def internals(fn, save=None):
    with open(fn) as handle:
        mol = load(handle, fmt="xyz")
    ic = InternalCoords(mol)
    B_mat = ic.B_matrix(mol)
    print(B_mat)
    if save:
        fmt = "% 1.4f"
        np.savetxt(save, B_mat, fmt=fmt)

# H2O2
h2o2_fn = "/scratch/programme/pysisyphus/xyz_files/h2o2_hf_321g_opt.xyz"
internals(h2o2_fn, "h2o2_berny.bmat")
