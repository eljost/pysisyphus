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

THIS_DIR = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))

def prepare_geometry():
    keywords = "HF 4-22GSP TightSCF"
    xyz_fn = "01_irc_sn2_fluour_transfer_optts.xyz"
    #blocks = "%pal nprocs 3 end"
    blocks = ""

    atoms, coords = parse_xyz_file(THIS_DIR / xyz_fn)
    coords *= ANG2BOHR
    geometry = Geometry(atoms, coords.flatten())

    geometry.set_calculator(ORCA(keywords, charge=-1, mult=1, blocks=blocks))

    hessian = geometry.hessian
    return geometry, THIS_DIR

def param_plot(irc, prefix):
    p1inds = (4, 0)
    p2inds = (5, 0)
    param_plot = ParamPlot(irc.all_coords, p1inds, p2inds)
    param_plot.plot()
    param_plot.save(THIS_DIR, prefix)
    param_plot.save_coords(THIS_DIR, prefix)


@pytest.mark.orca_irc
def test_irc_sn2_ftransfer_gs():
    geometry, THIS_DIR = prepare_geometry()
    prefix = "sn2_ftransfer_422gsp_gs"
    irc = GonzalesSchlegel(geometry, max_steps=5, step_length=0.1)
    irc.run()
    irc.write_trj(THIS_DIR, prefix)
    param_plot(irc, prefix)


@pytest.mark.orca_irc
def test_irc_sn2_ftransfer_gs03():
    # Fails
    geometry, THIS_DIR = prepare_geometry()
    prefix = "sn2_ftransfer_422gsp_gs03"
    irc = GonzalesSchlegel(geometry, max_steps=15, step_length=0.3)
    irc.run()
    irc.write_trj(THIS_DIR, prefix)
    param_plot(irc, prefix)


@pytest.mark.orca_irc
def test_irc_sn2_ftransfer_dvv():
    geometry, THIS_DIR = prepare_geometry()
    prefix = "sn2_ftransfer_422gsp_dvv"
    irc = DampedVelocityVerlet(geometry, max_steps=50)
    irc.run()
    irc.write_trj(THIS_DIR, prefix)
    param_plot(irc, prefix)


@pytest.mark.orca_irc
def test_irc_sn2_ftransfer_euler():
    geometry, THIS_DIR = prepare_geometry()
    prefix = "sn2_ftransfer_422gsp_euler"
    irc = Euler(geometry, max_steps=50)
    irc.run()
    irc.write_trj(THIS_DIR, prefix)
    param_plot(irc, prefix)

if __name__ == "__main__":
    test_irc_sn2_ftransfer_gs()
    #test_irc_sn2_ftransfer_gs03()
    #test_irc_sn2_ftransfer_dvv()
    #test_irc_sn2_ftransfer_euler()
    
