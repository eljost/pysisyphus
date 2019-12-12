#!/usr/bin/env python3

import os
import pathlib
from pathlib import Path

import numpy as np

from pysisyphus.helpers import geom_from_library, geom_from_xyz_file
from pysisyphus.calculators.Turbomole import Turbomole

THIS_DIR = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))


def test_forces():
    geom = geom_from_xyz_file(THIS_DIR / "h2o.xyz")
    control_path = THIS_DIR / "h2o_forces"
    turbo = Turbomole(control_path)
    geom.set_calculator(turbo)
    forces = geom.forces
    ref_forces = -np.loadtxt("h2o_gradient.ref")
    np.testing.assert_allclose(forces, ref_forces, atol=1e-8)


def test_benzene_wfoverlap():
    control_path = THIS_DIR / "benzene_3states"
    td_vec_fn = str(control_path / "ciss_a")
    mos_fn = str(control_path / "mos")
    turbo = Turbomole(control_path)
    geom = geom_from_library("benzene_bp86sto3g_opt.xyz")
    geom.set_calculator(turbo)
    turbo.td_vec_fn = td_vec_fn
    turbo.mos = mos_fn


def test_h2o_wfoverlap():
    path1 = THIS_DIR / "h2o1_wfoverlap"
    path2 = THIS_DIR / "h2o2_wfoverlap"
    geom1 = geom_from_xyz_file(path1 / "h2o1.xyz")
    geom2 = geom_from_xyz_file(path2 / "h2o2.xyz")
    turbo = Turbomole(path1)
    get_mos_ciss = lambda path: (str(path / "mos"), str(path / "ciss_a"))
    mos1, ciss1 = get_mos_ciss(path1)
    mos2, ciss2 = get_mos_ciss(path2)

    turbo.td_vec_fn = ciss1
    turbo.mos = mos1

    turbo.td_vec_fn = ciss2
    turbo.mos = mos2


if __name__ == "__main__":
    test_forces()
    #test_benzene_wfoverlap()
    #test_h2o_wfoverlap()
