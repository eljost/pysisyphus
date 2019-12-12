#!/usr/bin/env python3
import copy

import pytest

from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.calculators.IDPP import idpp_interpolate
from pysisyphus.cos.NEB import NEB
from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.FIRE import FIRE
from pysisyphus.optimizers.BFGS import BFGS

from qchelper.geometry import parse_xyz_file

KWARGS = {
    "images": 5,
    "max_cycles": 20,
    "convergence": {
        "max_force_thresh": 2e-3,
    },
    "dt": 0.5
}

def get_geoms(images=KWARGS["images"]):
    initial = "xyz_files/hcn.xyz"
    ts_guess = "xyz_files/hcn_iso_ts.xyz"
    final = "xyz_files/nhc.xyz"
    xyz_fns = (initial, ts_guess, final)
    atoms_coords = [parse_xyz_file(fn) for fn in xyz_fns]
    geoms = [Geometry(atoms, coords.flatten()) for atoms, coords in atoms_coords]
    geoms = idpp_interpolate(geoms, images_between=images)

    return geoms, copy.copy(KWARGS)

def run_cos_opt(cos, Opt, images, **kwargs):
    opt = Opt(cos, **kwargs)
    for img in cos.images:
        img.set_calculator(ORCA())
    opt.run()

    return opt

@pytest.mark.skip()
def test_fire():
    geoms, kwargs = get_geoms()
    neb = NEB(geoms)
    opt = run_cos_opt(neb, FIRE, **kwargs)

    #assert(opt.is_converged)

@pytest.mark.skip()
def test_bfgs():
    geoms, kwargs = get_geoms()
    neb = NEB(geoms)
    opt = run_cos_opt(neb, BFGS, **kwargs)

    #assert(opt.is_converged)

if __name__ == "__main__":
    #test_bfgs()
    test_fire()
