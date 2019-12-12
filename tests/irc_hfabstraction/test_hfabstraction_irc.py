#!/usr/bin/env python3
import logging
import os
import pathlib

import numpy as np
import pytest

from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.constants import ANG2BOHR
from pysisyphus.Geometry import Geometry
from pysisyphus.irc.GonzalesSchlegel import GonzalesSchlegel
from pysisyphus.irc.DampedVelocityVerlet import DampedVelocityVerlet
from pysisyphus.irc.Euler import Euler
from pysisyphus.irc.ParamPlot import ParamPlot

from qchelper.geometry import parse_xyz_file

def prepare_geometry(keywords=None, xyz_fn=None):
    this_dir = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))

    if not keywords:
        keywords = "HF STO-3G TightSCF"
    if not xyz_fn:
        xyz_fn = "hfabstraction_sto3g_ts.xyz"
    #blocks = "%pal nprocs 3 end"
    blocks = ""

    atoms, coords = parse_xyz_file(this_dir / xyz_fn)
    coords *= ANG2BOHR
    geometry = Geometry(atoms, coords.flatten())

    geometry.set_calculator(ORCA(keywords, charge=0, mult=1, blocks=blocks))

    hessian = geometry.hessian
    return geometry, this_dir


@pytest.mark.orca_irc
def test_hfabstraction_iso_gs_hfsto3g():
    geometry, this_dir = prepare_geometry()
    irc = GonzalesSchlegel(geometry, max_steps=5, step_length=0.3)
    irc.run()
    irc.write_trj(this_dir)


@pytest.mark.orca_irc
def test_hfabstraction_iso_dvv_hfsto3g():
    geometry, this_dir = prepare_geometry()
    irc = DampedVelocityVerlet(geometry, max_steps=5, v0=0.2)
    irc.run()
    irc.write_trj(this_dir)


@pytest.mark.orca_irc
def test_hfabstraction_iso_euler_hfsto3g():
    geometry, this_dir = prepare_geometry()
    irc = Euler(geometry, max_steps=100, step_length=0.025)
    #irc = Euler(geometry, max_steps=2, step_length=0.01, mass_weight=False)
    irc.run()
    irc.write_trj(this_dir, "hf_sto3g_mw")


@pytest.mark.orca_irc
def test_hfabstraction_iso_euler_hfsto3g_no_mw():
    geometry, this_dir = prepare_geometry()
    irc = Euler(geometry, max_steps=100, mass_weight=False, step_length=0.025)
    irc.run()
    irc.write_trj(this_dir, "hf_sto3g_nomw")


@pytest.mark.orca_irc
def test_hfabstraction_iso_euler_hf422gsp():
    xyz_fn = "07_hfabstraction_hf422gsp_ts.xyz"
    keywords = "HF 4-22GSP TightSCF"
    geometry, this_dir = prepare_geometry(keywords, xyz_fn)
    irc = Euler(geometry, max_steps=175, step_length=0.05)
    irc.run()
    irc.write_trj(this_dir, "hf_422gsp_mw")


@pytest.mark.orca_irc
def test_hfabstraction_iso_euler_hf422gsp_no_mw():
    xyz_fn = "07_hfabstraction_hf422gsp_ts.xyz"
    keywords = "HF 4-22GSP TightSCF"
    prefix = "hf_422gsp_nomw"
    geometry, this_dir = prepare_geometry(keywords, xyz_fn)
    irc = Euler(geometry, max_steps=125, mass_weight=False, step_length=0.025)
    irc.run()
    irc.write_trj(this_dir, prefix)

    p1inds = (3, 7)
    p2inds = (0, 3, 7)
    param_plot = ParamPlot(irc.all_coords, p1inds, p2inds)
    param_plot.plot()
    param_plot.show()
    param_plot.save(this_dir, prefix)


@pytest.mark.orca_irc
def tmp():
    xyz_fn = "07_hfabstraction_hf422gsp_ts.xyz"
    keywords = "HF 4-22GSP TightSCF"
    geometry, this_dir = prepare_geometry(keywords, xyz_fn)
    irc = Euler(geometry, max_steps=1, step_length=0.025, forward=True, mass_weight=False)
    irc.run()
    irc.write_trj(this_dir, "hf_422gsp_mw")


if __name__ == "__main__":
    #test_hfabstraction_iso_gs_hfsto3g()
    #test_hfabstraction_iso_dvv_hfsto3g()
    #test_hfabstraction_iso_euler_hfsto3g()
    #test_hfabstraction_iso_euler_hfsto3g_no_mw()
    test_hfabstraction_iso_euler_hf422gsp()
    test_hfabstraction_iso_euler_hf422gsp_no_mw()
    #tmp()
