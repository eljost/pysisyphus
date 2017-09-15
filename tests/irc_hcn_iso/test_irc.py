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
from pysisyphus.irc.ParamPlot import ParamPlot

from qchelper.geometry import parse_xyz_file

THIS_DIR = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))

def prepare_geometry(xyz_fn, keywords):

    #blocks = "%pal nprocs 3 end"
    blocks = ""

    atoms, coords = parse_xyz_file(THIS_DIR / xyz_fn)
    coords *= ANG2BOHR
    geometry = Geometry(atoms, coords.flatten())

    geometry.set_calculator(ORCA(keywords, charge=0, mult=1, blocks=blocks))

    hessian = geometry.hessian
    return geometry


def param_plot(irc, prefix):
    p1inds = (0, 1)
    p2inds = (1, 0, 2)
    param_plot = ParamPlot(irc.all_coords, p1inds, p2inds)
    param_plot.plot()
    param_plot.save(THIS_DIR, prefix)
    param_plot.save_coords(THIS_DIR, prefix)


@pytest.mark.orca_irc
def test_hcn_iso_gs_hfsto3g():
    keywords = "HF STO-3G TightSCF"
    xyz_fn = "hcn_gs_hfsto3g_ts.xyz"
    geometry = prepare_geometry(xyz_fn, keywords)
    prefix = "hfsto3g_gs"

    irc = GonzalesSchlegel(geometry, keywords, max_steps=7, step_length=0.15)
    irc.run()
    irc.write_trj(THIS_DIR, prefix)
    param_plot(irc, prefix)


@pytest.mark.orca_irc
def test_hcn_iso_dvv_hfsto3g():
    keywords = "HF STO-3G TightSCF"
    xyz_fn = "hcn_gs_hfsto3g_ts.xyz"
    geometry = prepare_geometry(xyz_fn, keywords)
    prefix = "hfsto3g_dvv"

    irc = DampedVelocityVerlet(geometry, max_steps=25)
    irc.run()
    irc.write_trj(THIS_DIR, prefix)
    param_plot(irc, prefix)


@pytest.mark.orca_irc
def test_hcn_iso_gs_bp86def2sv_p():
    """Crash in iter 6 forward."""
    keywords = "BP86 def2-SV(P) TightSCF"
    xyz_fn = "01_hcn_ts_bp86def2sv_p_ts.xyz"
    geometry = prepare_geometry(xyz_fn, keywords)
    prefix = "bp86def2sv_p_gs"

    irc = GonzalesSchlegel(geometry, keywords, step_length=0.15)
    irc.run()
    irc.write_trj(THIS_DIR, prefix)
    param_plot(irc, prefix)


@pytest.mark.orca_irc
def test_hcn_iso_gs_hfdef2sv_p():
    keywords = "HF def2-SV(P) TightSCF"
    xyz_fn = "02_hcn_ts_hfdef2sv_p_ts.xyz"
    geometry = prepare_geometry(xyz_fn, keywords)
    prefix = "hfdef2sv_p_gs"

    irc = GonzalesSchlegel(geometry, keywords, step_length=0.15)
    irc.run()
    irc.write_trj(THIS_DIR, prefix)
    param_plot(irc, prefix)


if __name__ == "__main__":
    #test_hcn_iso_gs_hfsto3g()
    #test_hcn_iso_dvv_hfsto3g()
    #test_hcn_iso_gs_bp86def2sv_p()
    test_hcn_iso_gs_hfdef2sv_p()
