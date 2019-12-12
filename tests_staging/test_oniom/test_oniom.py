#!/usr/bin/env python3

import numpy as np

from pysisyphus.calculators.ONIOM import ONIOM
from pysisyphus.helpers import geom_from_library, geom_from_xyz_file
from pysisyphus.optimizers.RFOptimizer import RFOptimizer


def test_acetaldehyd():
    calc_dict = {
        "high": {
            "type": "g16",
            "route": "b3lyp d95v",
            "pal": 4,
        },
        "low": {
            "type": "g16",
            "route": "hf sto-3g",
            "pal": 4,
        },
    }
    high_inds = (4,5,6)
    oniom = ONIOM(calc_dict, high_inds)

    geom = geom_from_library("acetaldehyd_oniom.xyz", coord_type="redund")
    geom.set_calculator(oniom)

    forces = geom.forces.reshape(-1, 3) # internal forces...
    forces = geom._forces.reshape(-1, 3)
    energy = geom.energy

    # print("energy")
    # print(f"{energy:.8f}")

    forces_str = np.array2string(forces, formatter={"float": lambda f: f"{f: .8f}",})
    # print("forces")
    # print(forces_str)

    from pysisyphus.optimizers.RFOptimizer import RFOptimizer
    rfo = RFOptimizer(geom, trust_max=.3, dump=True, thresh="gau")
    rfo.run()


def test_acetaldehyd_psi4_xtb():
    calc_dict = {
        "high": {
            "type": "pypsi4",
            "method": "scf",
            "basis": "sto-3g",
        },
        "low": {
            "type": "pyxtb",
        },
    }
    high_inds = (4,5,6)
    oniom = ONIOM(calc_dict, high_inds)

    geom = geom_from_library("acetaldehyd_oniom.xyz", coord_type="redund")
    geom.set_calculator(oniom)

    from pysisyphus.optimizers.RFOptimizer import RFOptimizer
    rfo = RFOptimizer(geom, trust_max=.3, dump=True, thresh="gau", line_search=True)
    rfo.run()


def test_biaryl_solvated():
    calc_dict = {
        "high": {
            "type": "g16",
            "route": "pm6",
            "pal": 4,
        },
        "low": {
            "type": "xtb",
            "pal": 4,
        },
    }
    high_inds = list(range(30))
    oniom = ONIOM(calc_dict, high_inds)

    geom = geom_from_xyz_file("bare_solvated.xyz")
    geom.set_calculator(oniom)

    opt_kwargs = {
        # "trust_max": .3,
        "dump": True,
        # "thresh": "gau",
        "prefix": "pm6_biaryl_",
        "max_cycles": 200,
    }
    rfo = RFOptimizer(geom, **opt_kwargs)
    rfo.run()


if __name__ == "__main__":
    # test_acetaldehyd()
    test_acetaldehyd_psi4_xtb()
    # test_biaryl_solvated()
