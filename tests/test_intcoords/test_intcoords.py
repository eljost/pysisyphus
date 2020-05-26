#!/usr/bin/env python3

import numpy as np
import pytest
from pytest import approx

from pysisyphus.helpers import geom_from_library, geom_loader
from pysisyphus.intcoords.fragments import get_fragments
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.testing import using_pyscf


def numhess(geom, step_size=0.0001):
    coords = geom.coords
    cnum = len(coords)
    H = list()
    for i in range(cnum):
        print(f"Step {i+1}/{cnum}")
        step = np.zeros_like(coords)
        step[i] += step_size

        # Plus
        pl_coords = coords + step
        geom.coords = pl_coords
        pl_forces = geom.forces

        # Minus
        min_coords = coords - step
        geom.coords = min_coords
        min_forces = geom.forces

        fd = -(pl_forces - min_forces) / (2*step_size)
        H.append(fd)
    H = np.array(H)
    # Symmetrize
    H = (H+H.T)/2
    return H


def compare_hessians(ref_H, num_H, ref_rms):
    print("Reference hessian")
    print(ref_H)
    print("Findiff hessian")
    print(num_H)

    rms = np.sqrt(np.mean((ref_H-num_H)**2))
    print(f"rms(diff)={rms:.8f}")

    return rms == approx(ref_rms, abs=1e-6)


@using_pyscf
@pytest.mark.parametrize(
    "xyz_fn, coord_type, ref_rms", [
        ("hcn_bent.xyz", "cart", 1.2e-6),
        ("h2o2_rot2.xyz", "redund", 0.000876695),
    ]
)
def test_numhess(xyz_fn, coord_type, ref_rms):
    # geom = geom_from_library("hcn_bent.xyz")
    geom = geom_from_library(xyz_fn, coord_type=coord_type)

    # Interestingly enough the test will fail with keep_chk == True ...
    # as the RMS value will be much higher
    calc = PySCF(basis="321g", pal=2, keep_chk=False)
    geom.set_calculator(calc)

    H = geom.hessian
    nH = numhess(geom)

    assert compare_hessians(H, nH, ref_rms)


def test_get_fragments():
    geom = geom_loader("lib:thr75_from_1bl8.xyz")

    fragments = get_fragments(geom.atoms, geom.coords3d)
    assert len(fragments) == 4
