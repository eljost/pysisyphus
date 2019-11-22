#!/usr/bin/env python3

import numpy as np

from pysisyphus.init_logging import init_logging
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.calculators.Turbomole import Turbomole
from pysisyphus.optimizers.RFOptimizer import RFOptimizer
from pysisyphus.helpers import geom_from_library, geom_from_xyz_file


init_logging()


np.set_printoptions(suppress=True, precision=4)

def test_pyscf():
    # pyscf_ = PySCF(basis="3-21g", pal=4)
    pyscf_ = PySCF(basis="3-21g", pal=4)
    # pyscf_ = PySCF(basis="3-21g", method="mp2", pal=4)
    # geom = geom_from_library("birkholz/vitamin_c.xyz")
    # geom = geom_from_library("hcn.xyz")
    geom = geom_from_library("hcn_iso_ts.xyz")
    geom.set_calculator(pyscf_)
    f = geom.forces.reshape(-1, 3)
    print("PySCF")
    print(f)

    H = geom.hessian
    print("PySCF hessian")
    print(H.reshape(-1, 9))

    ref_geom = geom.copy()
    from pysisyphus.calculators.Gaussian16 import Gaussian16
    g16 = Gaussian16("hf/3-21G", pal=4)
    # g16 = Gaussian16("mp2/3-21G", pal=4)
    ref_geom.set_calculator(g16)
    f_ref = ref_geom.forces.reshape(-1, 3)
    print("Gaussian16")
    print(f_ref)

    # np.testing.assert_allclose(f, f_ref, rtol=5e-3)

    H_ref = ref_geom.hessian
    print("G16 Hess")
    print(H_ref)


def test_pyscf_tddft():
    calc_kwargs = {
        "xc": "b3lyp",
        "method": "tddft",
        "basis": "3-21g",
        "pal": 4,
        "nstates": 10,
        "root": 2,
        "track": True,
    }
    pyscf_ = PySCF(**calc_kwargs)
    geom = geom_from_library("hcn_iso_ts.xyz")
    geom.set_calculator(pyscf_)
    e = geom.energy
    print(e)
    f = geom.forces.reshape(-1, 3)
    print("PySCF")
    print(f)


def test_pyscf_cytosin_td_opt():
    calc_kwargs = {
        "xc": "pbe0",
        "method": "tddft",
        # "method": "dft",
        "basis": "def2SVP",
        # "auxbasis": "weigend",
        # "basis": "sto3g",
        "pal": 4,
        "nstates": 2,
        "root": 1,
        "track": True,
        "ovlp_with": "adapt",
        "ovlp_type": "tden",
    }
    pyscf_ = PySCF(**calc_kwargs)
    geom = geom_from_library("cytosin.xyz", coord_type="redund")
    geom.set_calculator(pyscf_)
    opt = RFOptimizer(geom)
    opt.run()


def turbo_comp():
    calc_kwargs = {
        "xc": "pbe",
        "method": "tddft",
        "basis": "def2-SVP",
        "pal": 4,
        "nstates": 2,
        "root": 1,
        "track": True,
    }
    pyscf_ = PySCF(**calc_kwargs)
    geom = geom_from_xyz_file("cytosin.xyz")
    geom.set_calculator(pyscf_)
    f = geom.forces.reshape(-1, 3)
    print("pyscf")
    print(f.reshape(-1,3))

    ref_geom = geom.copy()
    turbo_kwargs = {
            "control_path": "/scratch/turbontos/control_pbe_clean",
            "track": True,
            "ovlp_type": "wf",
            "ovlp_with": "first",
            "pal": 4,
    }
    turbo = Turbomole(**turbo_kwargs)
    ref_geom.set_calculator(turbo)
    f_ref = ref_geom.forces
    print("Turbo")
    print(f_ref.reshape(-1,3))


if __name__ == "__main__":
    # test_pyscf()
    # test_pyscf_tddft()
    test_pyscf_cytosin_td_opt()
    # turbo_comp()
