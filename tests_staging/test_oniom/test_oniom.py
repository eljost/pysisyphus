#!/usr/bin/env python3

import numpy as np
import pytest

from pysisyphus.calculators.ONIOM import ONIOM
from pysisyphus.calculators.ONIOMext import ONIOMext
from pysisyphus.helpers import geom_from_library, geom_from_xyz_file
from pysisyphus.optimizers.RFOptimizer import RFOptimizer

from pysisyphus.init_logging import init_logging


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

    assert geom.energy == pytest.approx(-153.07432042299052)
    assert np.linalg.norm(geom.forces) == pytest.approx(0.03768246934785125)

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


def test_oniomext():
    geom = geom_from_library("alkyl17_sto3g_opt.xyz")

    real = set(range(len(geom.atoms)))
    medmin = set((0,1,2,3,4,5,6, 46,47,48,49,50,51,52))
    med = list(real - medmin)
    h1 = list(range(13, 22))
    h2 = list(range(31, 40))

    calcs = {
        "real": {
            "route": "HF/STO-3G",
        },
        "medium": {
            "route": "HF/3-21G",
        },
        "high1": {
            "route": "HF/6-31G",
        },
        "high2": {
            "route": "HF/6-311G",
        },
    }
    for key, calc in calcs.items():
        calc["type"] = "g16"
        calc["pal"] = 2
        calc["mult"] = 1
        calc["charge"] = 0

    models = {
        "med" : {
            # "inds": (2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
            "inds": med,
            "calc": "medium",
        },
        "h1": {
            # "inds": (4, 5, 6),
            "inds": h1,
            "calc": "high1",
        },
        "h2": {
            # "inds": (10, 11, 12),
            "inds": h2,
            "calc": "high2",
        }
    }

    layers = [["h1", "h2"], "med"]

    oniom = ONIOMext(calcs, models, geom, layers)

    assert oniom.layer_num == 3

    geom.set_calculator(oniom)

    en = geom.energy
    print("energy", en)
    assert en == pytest.approx(-661.3512410069466)


def test_oniomext_gradient():
    geom = geom_from_library("acetaldehyd_oniom.xyz", coord_type="redund")

    real = set(range(len(geom.atoms)))
    high = [4, 5, 6]

    calcs = {
        "real": {
            "route": "hf sto-3g",
        },
        "high": {
            "route": "b3lyp d95v",
        },
    }
    for key, calc in calcs.items():
        calc["type"] = "g16"
        calc["pal"] = 4
        calc["mult"] = 1
        calc["charge"] = 0

    models = {
        "high": {
            "inds": high,
            "calc": "high",
        },
    }

    # layers = ["high"]
    layers = None

    oniom = ONIOMext(calcs, models, geom, layers)

    assert oniom.layer_num == 2

    geom.set_calculator(oniom)

    # Calculate forces
    assert np.linalg.norm(geom.forces) == pytest.approx(0.03768246934785125)
    assert geom.energy == pytest.approx(-153.07432042299052)


def test_oniomext_ee():
    geom = geom_from_library("oniom_ee_model_system.xyz", coord_type="redund")

    all_ = set(range(len(geom.atoms)))
    high = list(sorted(all_ - set((21, 20, 19, 15, 14, 13))))

    calcs = {
        "real": {
            "route": "hf 6-31G",
        },
        "high": {
            "route": "mp2 6-31G*",
        },
    }
    for key, calc in calcs.items():
        calc["type"] = "g16"
        calc["pal"] = 4
        calc["mult"] = 1
        calc["charge"] = 0

    models = {
        "high": {
            "inds": high,
            "calc": "high",
        },
    }

    oniom = ONIOMext(calcs, models, geom)
    geom.set_calculator(oniom)

    # Calculate forces and energy
    forces = geom.forces
    energy = geom.energy

    assert np.linalg.norm(forces) == pytest.approx(0.0585419180)
    assert energy == pytest.approx(-588.02530947)


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
    init_logging()
    # test_acetaldehyd()
    # test_oniomext()
    # test_oniomext_gradient()
    test_oniomext_ee()
    # test_acetaldehyd_psi4_xtb()
    # test_biaryl_solvated()
