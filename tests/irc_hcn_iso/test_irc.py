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
#grad = -geometry.forces
#energy = geometry.energy
#print("energy", energy)
#print("gradient", grad)
hessian = geometry.hessian
#print("hessian", hessian)
#np.save(grad_fn, grad)
#np.save(hessian_fn, hessian)

#gradient = np.load(grad_fn)
#geometry.gradient = gradient
#hessian = np.load(hessian_fn)
#geometry.hessian = hessian

gs_irc = GonzalesSchlegel(geometry, max_cycles=10)
gs_irc.run()
xyz_strings = list()
for irc_coords in gs_irc.coords:
    geometry.coords = irc_coords
    as_xyz = geometry.as_xyz()
    xyz_strings.append(as_xyz)

xyzs_joined = "\n".join(xyz_strings)
with open(this_dir / "hcn_iso_irc.trj", "w") as handle:
    handle.write(xyzs_joined)

np.savetxt(this_dir / "hcn_iso_irc.energies", gs_irc.energies)

#coords_save = np.savetxt("hcn_iso_coords.dat", gs_irc.coords)
