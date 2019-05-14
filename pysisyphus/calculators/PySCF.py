#!/usr/bin/env python3

import pyscf
from pyscf import gto, grad, lib

from pysisyphus.calculators.OverlapCalculator import OverlapCalculator


class PySCF(OverlapCalculator):

    conf_key = "pyscf"
    drivers = {
        ("dft", False): pyscf.dft.RKS,
        ("dft", True): pyscf.dft.UKS,
        ("scf", False): pyscf.scf.RHF,
        ("scf", True): pyscf.scf.UHF,
    }
    # grad_drivers = {
        # ("dft", False): pyscf.grad.RKS,
        # ("dft", True): pyscf.grad.UKS,
        # ("scf", False): pyscf.grad.RHF,
        # ("scf", True): pyscf.grad.UHF,
    # }

    def __init__(self, method, basis, xc=None, mem=2000, **kwargs):
        super().__init__(**kwargs)

        self.method = method
        self.basis = basis
        self.xc = xc
        self.mem = mem

        self.out_fn = "pyscf.out"
        self.unrestricted = self.mult > 1
        # self.to_keep = ("pyscf.inp", "pyscf.out", )

        lib.num_threads(self.pal)

    def get_driver(self, mol):
        if self.xc:
            driver = self.drivers[("dft", self.unrestricted)]
            driver.xc = self.xc
            mf = driver(mol)
        else:
            driver = self.drivers[("scf", self.unrestricted)]
            mf = driver(mol)
        return mf

    # def get_grad_driver(self):
        # if self.xc:
            # driver = self.grad_drivers[("dft", self.unrestricted)]
        # else:
            # driver = self.grad_drivers[("scf", self.unrestricted)]
        # return driver

    def prepare_input(self, atoms, coords):
        mol = gto.Mole()
        mol.atom = [(atom, c) for atom, c in zip(atoms, coords)]
        mol.basis = self.basis
        mol.unit = "Bohr"
        mol.charge = self.charge
        mol.spin = self.mult - 1
        mol.symmetry = False
        mol.verbose = 5
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
    from pysisyphus.helpers import geom_from_library
    from pysisyphus.Geometry import Geometry
    from pysisyphus.optimizers.RFOptimizer import RFOptimizer

    pyscf_ = PySCF(method="dft", xc="b3lyp", basis="631g", pal=8)
    # geom = geom_from_library("birkholz/artemisin.xyz")
    geom = geom_from_library("hcn.xyz")
    print(geom)
    geom.set_calculator(pyscf_)
    geom.forces
    opt = RFOptimizer(geom)
    opt.run()
