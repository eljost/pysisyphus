#!/usr/bin/env python3
import logging
import os
import pathlib

import numpy as np
import pytest

from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.Geometry import Geometry
from pysisyphus.irc.GonzalesSchlegel import GonzalesSchlegel

from qchelper.geometry import parse_xyz_file

np.set_printoptions(suppress=True, precision=4)

#https://verahill.blogspot.de/2013/06/439-calculate-frequencies-from-hessian.html
#https://chemistry.stackexchange.com/questions/74639

this_dir = pathlib.Path(os.path.dirname(os.path.realpath(__file__)))
grad_fn = this_dir / "hcn_gradient.npy"
hessian_fn = this_dir / "hcn_hessian.npy"
initial = this_dir / "hcn_ts_standard.xyz"

atoms, coords = parse_xyz_file(initial)
geometry = Geometry(atoms, coords.flatten())
geometry.set_calculator(ORCA())
#grad = -geom.forces
#energy = geom.energy
#print("energy", energy)
#print("gradient", grad)
#hessian = geom.hessian
#print("hessian", hessian)
#np.save(grad_fn, grad)
#np.save(hessian_fn, hessian)

gradient = np.load(grad_fn)
geometry.gradient = gradient
hessian = np.load(hessian_fn)
geometry.hessian = hessian

gs_irc = GonzalesSchlegel(geometry)
