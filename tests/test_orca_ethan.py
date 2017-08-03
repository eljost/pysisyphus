#!/usr/bin/env python3

import copy

import pytest

from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.Geometry import Geometry
from pysisyphus.optimizers.FIRE import FIRE
from pysisyphus.optimizers.BFGS import BFGS

from qchelper.geometry import parse_xyz_file

KWARGS = {
    "max_cycles": 30,
    "max_force_thresh": 0.002,
    #"dt": 0.5,
}

def prepare_opt():
    ethan_xyz = "xyz_files/ethan.xyz"
    atoms, coords = parse_xyz_file(ethan_xyz)
    geom = Geometry(atoms, coords.flatten())
    geom.set_calculator(ORCA())
    
    return geom, copy.copy(KWARGS)

def run_opt(geom, Optimizer, **kwargs):
    opt = Optimizer(geom, **kwargs)
    opt.run()
    
    return opt

@pytest.mark.skip()
def test_fire():
    geom, kwargs = prepare_opt()
    opt = run_opt(geom, FIRE, **kwargs)

    assert(opt.is_converged)

@pytest.mark.skip()
def test_bfgs():
    geom, kwargs = prepare_opt()
    kwargs["max_force_thresh"] = 0.0005
    opt = run_opt(geom, BFGS, **kwargs)
    #print(opt.energies[-1]*27.211386)

    assert(opt.is_converged)
 
if __name__ == "__main__":
    test_fire()
    test_bfgs()
