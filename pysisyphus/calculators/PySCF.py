#!/usr/bin/env python3

import pyscf
from pyscf import gto, grad, lib

from pysisyphus.calculators.OverlapCalculator import OverlapCalculator


class PySCF(OverlapCalculator):

    conf_key = "pyscf"
    drivers = {
        # key: (method, unrestricted?)
        ("dft", False): pyscf.dft.RKS,
        ("dft", True): pyscf.dft.UKS,
        ("scf", False): pyscf.scf.RHF,
        ("scf", True): pyscf.scf.UHF,
        # ("mp2", False): pyscf.mp.MP2,
        # ("mp2", True): pyscf.mp.UMP2,
    }

    def __init__(self, basis, xc=None, method=None,  mem=2000, **kwargs):
        super().__init__(**kwargs)

        self.basis = basis
        self.xc = xc
        self.method = method
        self.mem = mem

        self.out_fn = "pyscf.out"
        self.unrestricted = self.mult > 1

        lib.num_threads(self.pal)

    def get_driver(self, mol):
        if self.xc:
            driver = self.drivers[("dft", self.unrestricted)]
            driver.xc = self.xc
            mf = driver(mol)
        # elif self.method == "mp2":
            # hf_driver = self.drivers[("scf", self.unrestricted)]
            # mf = hf_driver(mol)
            # mp2_driver = self.drivers[("mp2", self.unrestricted)]
            # mf = mp2_driver(mf)
        else:
            driver = self.drivers[("scf", self.unrestricted)]
            mf = driver(mol)
        mf.conv_tol = 1e-8
        mf.max_cycle = 150
        return mf

    def prepare_input(self, atoms, coords):
        mol = gto.Mole()
        mol.atom = [(atom, c) for atom, c
                    in zip(atoms, coords.reshape(-1, 3))]
        mol.basis = self.basis
        mol.unit = "Bohr"
        mol.charge = self.charge
        mol.spin = self.mult - 1
        mol.symmetry = False
        mol.verbose = 4
        mol.output = self.out_fn
        mol.max_memory = self.mem * self.pal
        mol.build()

        return mol

    def get_energy(self, atoms, coords):
        mol = self.prepare_input(atoms, coords)
        mf = self.get_driver(mol)
        mf.kernel()
        results = {
            "energy": mf.energy_tot(),
        }
        return results

    def get_forces(self, atoms, coords):
        mol = self.prepare_input(atoms, coords)
        mf = self.get_driver(mol)
        # >>> mf.chkfile = '/path/to/chkfile'
        # >>> mf.init_guess = 'chkfile'
        mf.kernel()
        grad_driver = mf.Gradients()
        gradient = grad_driver.kernel()

        results = {
            "energy": mf.energy_tot(),
            "forces": -gradient.flatten(),
        }

        return results

    def __str__(self):
        return f"PySCF({self.name})"


if __name__ == "__main__":
    import numpy as np
    np.set_printoptions(suppress=True, precision=4)
    from pysisyphus.helpers import geom_from_library

    # pyscf_ = PySCF(method="mp2", basis="3-21g", pal=4)
    pyscf_ = PySCF(basis="3-21g", pal=4)
    geom = geom_from_library("birkholz/vitamin_c.xyz")
    geom.set_calculator(pyscf_)
    f = geom.forces.reshape(-1, 3)
    print("PySCF")
    print(f)

    ref_geom = geom.copy()
    from pysisyphus.calculators.Gaussian16 import Gaussian16
    g16 = Gaussian16("hf/3-21G", pal=4)
    # g16 = Gaussian16("mp2/3-21G", pal=4)
    ref_geom.set_calculator(g16)
    f_ref = ref_geom.forces.reshape(-1, 3)
    print("Gaussian16")
    print(f_ref)

    np.testing.assert_allclose(f, f_ref, rtol=5e-3)
