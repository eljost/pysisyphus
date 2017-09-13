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

from qchelper.geometry import parse_xyz_file

def run_irc(xyz_fn, keywords):
    this_dir = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))

    #blocks = "%pal nprocs 3 end"
    blocks = ""

    atoms, coords = parse_xyz_file(this_dir / xyz_fn)
    coords *= ANG2BOHR
    geometry = Geometry(atoms, coords.flatten())

    geometry.set_calculator(ORCA(keywords, charge=0, mult=1, blocks=blocks))

    hessian = geometry.hessian
    gs_irc = GonzalesSchlegel(geometry, max_steps=2, step_length=0.1, backward=True)

    gs_irc.run()

    return gs_irc

@pytest.mark.orca_irc
def test_hcn_iso_hfsto3g():
    keywords = "HF STO-3G TightSCF"
    xyz_fn = "hcn_gs_hfsto3g_ts.xyz"
    irc = run_irc(xyz_fn, keywords)


@pytest.mark.orca_irc
def test_hcn_iso_bp86def2sv_p():
    keywords = "BP86 def2-SV(P) TightSCF"
    xyz_fn = "01_hcn_ts_bp86def2sv_p_ts.xyz"
    irc = run_irc(xyz_fn, keywords)


@pytest.mark.orca_irc
def test_hcn_iso_hfdef2sv_p():
    keywords = "HF def2-SV(P) TightSCF"
    xyz_fn = "02_hcn_ts_hfdef2sv_p_ts.xyz"
    irc = run_irc(xyz_fn, keywords)


if __name__ == "__main__":
    test_hcn_iso_hfsto3g()
    #test_hcn_iso_bp86def2sv_p()
    #test_hcn_iso_hfdef2sv_p()
