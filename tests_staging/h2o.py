#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from calculators.ORCA import ORCA
from cos.NEB import NEB
from Geometry import Geometry
from optimizers.SteepestDescent import SteepestDescent
from optimizers.FIRE import FIRE

from qchelper.geometry import parse_xyz_file

CYCLES = 15
IMAGES = 3

def get_geoms():
    educt = "xyz_files/h2o_inv_educt.xyz"
    product = "xyz_files/h2o_inv_product.xyz" 
    xyz_fns = (educt, product)
    atoms_coords = [parse_xyz_file(fn) for fn in xyz_fns]
    geoms = [Geometry(atoms, coords.flatten()) for atoms, coords in atoms_coords]
    return geoms


def run_cos_opt(cos):
    cos.interpolate(IMAGES)
    for img in cos.images:
        img.set_calculator(ORCA())

    kwargs = {
        #"rel_step_thresh": 1e-3,
        "max_cycles": CYCLES,
    }
    opt = FIRE(cos, **kwargs)
    #opt = SteepestDescent(cos, **kwargs)
    opt.run()


if __name__ == "__main__":
    """
    try:
        neb = NEB(get_geoms())
        run_cos_opt(neb)
    except:
        import pdb
        pdb.set_trace()
    """
    neb = NEB(get_geoms())
    run_cos_opt(neb)
