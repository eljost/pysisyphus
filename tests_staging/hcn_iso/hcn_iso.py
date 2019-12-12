#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np

from calculators.ORCA import ORCA
from calculators.IDPP import idpp_interpolate
from cos.NEB import NEB
from cos.SimpleZTS import SimpleZTS
from Geometry import Geometry
from optimizers.SteepestDescent import SteepestDescent
from optimizers.FIRE import FIRE

from qchelper.geometry import parse_xyz_file

CYCLES = 5
IMAGES = 3


def get_geoms():
    educt = "xyz_files/hcn.xyz"
    ts_guess ="xyz_files/hcn_iso_ts.xyz"
    product = "xyz_files/nhc.xyz" 
    xyz_fns = (educt, ts_guess, product)
    #xyz_fns = (educt, ts_guess)
    atoms_coords = [parse_xyz_file(fn) for fn in xyz_fns]
    geoms = [Geometry(atoms, coords.flatten()) for atoms, coords in atoms_coords]
    return geoms


def get_idpp_geoms():
    geoms = get_geoms()
    return idpp_interpolate(geoms, images_between=IMAGES)


def run_cos_opt(cos, idpp=False):
    # Otherwise it's already interpolated
    if not idpp:
        cos.interpolate(IMAGES)
    for img in cos.images:
        img.set_calculator(ORCA())

    kwargs = {
        "max_cycles":CYCLES,
    }
    opt = FIRE(cos, **kwargs)
    #opt = SteepestDescent(cos, **kwargs)
    opt.run()


def procrustes_test(cos):
    cos.interpolate(IMAGES)
    for img in cos.images:
        img.set_calculator(ORCA())

    kwargs = {
        "max_cycles":CYCLES,
    }
    opt = FIRE(cos, **kwargs)
    opt.procrustes()


if __name__ == "__main__":
    """
    try:
        neb = NEB(get_geoms())
        run_cos_opt(neb)
    except:
        import pdb
        pdb.set_trace()
    neb = NEB(get_geoms())
    run_cos_opt(neb)

    szts_equal = SimpleZTS(get_geoms(), param="equal")
    run_cos_opt(szts_equal)
    print()
    """
    #procrustes_test(neb)
    """
    # IDPP
    neb = NEB(get_idpp_geoms())
    run_cos_opt(neb, idpp=True)
    """
    get_idpp_geoms()

    """
    # Linear
    neb = NEB(get_geoms())
    run_cos_opt(neb)
    """
