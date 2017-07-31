#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from calculators.ORCA import ORCA
from cos.NEB import NEB
from Geometry import Geometry
from optimizers.SteepestDescent import SteepestDescent
from optimizers.FIRE import FIRE

from qchelper.geometry import parse_xyz_file

CYCLES = 5
IMAGES = 1

"""
def run_new():
    neb.save("iter000.trj")
    ics = list()
    for i in range(1, 7):
        neb.cycle()
        inner_coords = neb.images[1].coords
        ics.append(inner_coords)
        neb.save("iter{:03}.trj".format(i))

    ic_reshaped = [ic.reshape((-1,3)) for ic in ics]
    atoms = atoms_coords[0][0]
    from qchelper.geometry import make_trj_str
    trj_str = make_trj_str(atoms, ic_reshaped)
    with open("inner_coords.trj", "w") as handle:
        handle.write(trj_str)
"""

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

    #opt = SteepestDescent(cos,
    opt = FIRE(cos,
               max_cycles=CYCLES)
    opt.run()


if __name__ == "__main__":
    neb = NEB(get_geoms())
    run_cos_opt(neb)
