#!/usr/bin/env python3

from time import time

from pysisyphus.helpers import geom_from_library
from pysisyphus.cos.NEB import NEB
from pysisyphus.interpolate import interpolate_all
from pysisyphus.calculators import Psi4
from pysisyphus.optimizers.QuickMin import QuickMin
from pysisyphus.optimizers.FIRE import FIRE
from pysisyphus.optimizers.LBFGS import LBFGS


def test_psi4():
    geom = geom_from_library("hcn_iso_ts.xyz")
    psi4_kwargs = {
        "pal": 4,
        "mem": 2000,
        # "method": "hf",
        "method": "b3lyp",
        "basis": "def2-svp",
        # "mult": 1,
        # "charge": 2,
    }
    psi4 = Psi4(**psi4_kwargs)
    geom.set_calculator(psi4)
    print(psi4.base_cmd)
    # en = geom.energy
    # print(en)
    f = geom.forces
    print(f)
    e = geom.energy
    print(e)

    start = time()
    h = geom.hessian
    end = time()
    print(h)
    dur = end - start
    print("hess calc took", int(dur), "seconds")


def test_psi4_neb():
    start = geom_from_library("hcn.xyz")
    ts = geom_from_library("hcn_iso_ts_guess.xyz")
    end = geom_from_library("nhc.xyz")
    images = interpolate_all((start, ts, end), 4, kind="lst")
    for i, img in enumerate(images):
        img.set_calculator(Psi4("tpss", "sto-3g", pal=4, calc_number=i))
    neb = NEB(images, climb=True)
    with open("interpolated.trj", "w") as handle:
        handle.write(neb.as_xyz())
    # qm = QuickMin(neb, align=True, dump=True, max_cycles=10)
    # qm.run()
    fire = FIRE(neb, align=True, dump=True, max_cycles=10)
    fire.run()

    lbfgs = LBFGS(neb, align=True, dump=True)
    lbfgs.run()


if __name__ == "__main__":
    # test_psi4()
    test_psi4_neb()
