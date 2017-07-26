#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from calculators.ORCA import ORCA
from cos.NEB import NEB
from Geometry import Geometry
from optimizers.SteepestDescent import SteepestDescent

from qchelper.geometry import parse_xyz_file


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

def run_neb():
    educt = "h2o_inv_educt.xyz"
    product = "h2o_inv_product.xyz" 
    xyz_fns = [educt, product]
    atoms_coords = [parse_xyz_file(fn) for fn in xyz_fns]
    geoms = [Geometry(atoms, coords.flatten()) for atoms, coords in atoms_coords]
    neb = NEB(geoms)
    neb.interpolate_images(2)
    for img in neb.images[1:-1]:
        img.set_calculator(ORCA())

    sd = SteepestDescent(neb, max_cycles=5)
    sd.run()


def run_opt():
    atoms, coords = parse_xyz_file("h2o_inv_educt.xyz")
    geom = Geometry(atoms, coords.flatten())
    geom.set_calculator(ORCA())
    sd = SteepestDescent(geom, max_cycles=10)
    sd.run()

if __name__ == "__main__":
    run_neb()
    run_opt()

