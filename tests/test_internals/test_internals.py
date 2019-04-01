#!/usr/bin/env python3

from pysisyphus.constants import ANG2BOHR, BOHR2ANG
from pysisyphus.Geometry import Geometry
from pysisyphus.helpers import geom_from_library
from pysisyphus.calculators.XTB import XTB
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.optimizers.SteepestDescent import SteepestDescent

import numpy as np

import sys
sys.path.insert(0, "/home/carpx/Code/pyberny")
from berny import Berny, geomlib, optimize


def run():
    fn = "codein.xyz"
    # fn = "h2o2_guess.xyz"
    geom = geom_from_library(fn, coord_type="redund")
    # geom = geom_from_library(fn)#, coord_type="redund")
    # ints = geom.internal
    # bonds, bends, tors = ints.prim_indices
    # geom = geom_from_xyz_file(fn)
    xtb = XTB(acc=0.1)
    geom.set_calculator(xtb)
    opt_kwargs = {
        # "max_cycles": 3,
        "dump": True,
        # "thresh": "gau",
    }
    opt = RFOptimizer(geom, **opt_kwargs)
    # opt = SteepestDescent(geom, **opt_kwargs)
    opt.run()
    with open("final_pysis.xyz", "w") as handle:
        handle.write(geom.as_xyz())


def test_pyberny():
    fn = "codein.xyz"
    # fn = "h2o2_guess.xyz"

    def solver():
        atoms, _ = yield

        atoms, coords = zip(*atoms)
        coords = np.array(coords).flatten() * ANG2BOHR
        geom = Geometry(atoms, coords)
        xtb = XTB(acc=0.1)
        geom.set_calculator(xtb)

        while True:
            gradient = geom.gradient
            energy = geom.energy
            atoms, _ = yield energy, gradient
            _, coords = zip(*atoms)
            coords = np.array(coords).flatten() * ANG2BOHR
            geom.coords = coords

    berny = Berny(geomlib.readfile(fn))
    final = optimize(berny, solver())
    final.write("final_pyberny.xyz")


def test_pyberny_mod():
    pass


if __name__ == "__main__":
    run()
    # test_pyberny()
    # test_pyberny_mod()

