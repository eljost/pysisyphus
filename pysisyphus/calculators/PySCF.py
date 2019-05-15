#!/usr/bin/env python3

import numpy as np
import pyscf
from pyscf import gto, grad, lib, hessian

from pysisyphus.calculators.OverlapCalculator import OverlapCalculator


class PySCF(OverlapCalculator):

    conf_key = "pyscf"
    drivers = {
        # key: (method, unrestricted?)
        ("dft", False): pyscf.dft.RKS,
        ("dft", True): pyscf.dft.UKS,
        ("scf", False): pyscf.scf.RHF,
        ("scf", True): pyscf.scf.UHF,
        ("mp2", False): pyscf.mp.MP2,
        ("mp2", True): pyscf.mp.UMP2,
    }
    multisteps = {
        "scf": ("scf", ),
        "dft": ("dft", ),
        "mp2": ("scf", "mp2"),
        "tddft": ("dft", "tddft"),
    }

    def __init__(self, basis, xc=None, method="scf",  mem=2000, **kwargs):
        super().__init__(**kwargs)

        self.basis = basis
        self.xc = xc
        self.method = method
        if self.xc and self.method != "tddft":
            self.method = "dft"
        # self.multistep = self.method in ("mp2 tddft".split())
        self.mem = mem

        self.out_fn = "pyscf.out"
        self.unrestricted = self.mult > 1

        lib.num_threads(self.pal)

    def get_driver(self, step, mol=None, mf=None):
        if mol and (step == "dft"):
            driver = self.drivers[(step, self.unrestricted)]
            driver.xc = self.xc
            mf = driver(mol)
            mf.conv_tol = 1e-8
            mf.max_cycle = 150
        elif mol and (step == "scf"):
            driver = self.drivers[(step, self.unrestricted)]
            mf = driver(mol)
            mf.conv_tol = 1e-8
            mf.max_cycle = 150
        elif mf and (step == "mp2"):
            mp2_mf = self.drivers[(step, self.unrestricted)]
            mf = mp2_mf(mf)
        else:
            raise Exception("Unknown method!")
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
        mol.verbose = 5
        mol.output = self.out_fn
        mol.max_memory = self.mem * self.pal
        mol.build()

        return mol

    def get_energy(self, atoms, coords):
        mol = self.prepare_input(atoms, coords)
        mf = self.run(mol)
        results = {
            "energy": mf.e_tot,
        }
        return results

    def get_forces(self, atoms, coords):
        mol = self.prepare_input(atoms, coords)
        mf = self.run(mol)
        # >>> mf.chkfile = '/path/to/chkfile'
        # >>> mf.init_guess = 'chkfile'
        grad_driver = mf.Gradients()
        gradient = grad_driver.kernel()

        results = {
            "energy": mf.e_tot,
            "forces": -gradient.flatten(),
        }

        return results

    def get_hessian(self, atoms, coords):
        mol = self.prepare_input(atoms, coords)
        mf = self.run(mol)
        H = mf.Hessian().kernel()

        # The returned hessian is 4d ... ok. This probably serves a purpose
        # that I don't understand. We transform H to a nice 2d array.
        H = np.hstack(np.concatenate(H, axis=1))
        results = {
            "energy": mf.e_tot,
            "hessian": H,
        }
        return results

    def run(self, mol):
        steps = self.multisteps[self.method]
        self.log(f"Running steps '{steps}' for method {self.method}")
        for i, step in enumerate(steps):
            if i == 0:
                mf = self.get_driver(step, mol=mol)
            else:
                mf = self.get_driver(step, mf=prev_mf)
            mf.kernel()
            self.log(f"Completed {step} step")
            prev_mf = mf
        return mf

    def __str__(self):
        return f"PySCF({self.name})"
