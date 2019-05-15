#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pyscf
from pyscf import gto, grad, lib, hessian, tddft
from pyscf.dft import xcfun

from pysisyphus.calculators.OverlapCalculator import OverlapCalculator
from pysisyphus.calculators.WFOWrapper import WFOWrapper


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
        "tda": ("dft", "tda"),
    }

    def __init__(self, basis, xc=None, method="scf",  mem=2000,
                 root=None, nstates=None, **kwargs):
        super().__init__(**kwargs)

        self.basis = basis
        self.xc = xc
        self.method = method
        if self.xc and self.method != "tddft":
            self.method = "dft"
        # self.multistep = self.method in ("mp2 tddft".split())
        self.mem = mem
        self.root = root
        self.nstates = nstates
        if self.method == "tddft":
            assert self.nstates, "nstates must be set with method='tddft'!"
        if self.track:
            assert self.root <= self.nstates, "Please supply 'root' with " \
                "'track: True'!"

        self.out_fn = "pyscf.out"
        self.unrestricted = self.mult > 1

        lib.num_threads(self.pal)

    def get_driver(self, step, mol=None, mf=None):
        if mol and (step == "dft"):
            driver = self.drivers[(step, self.unrestricted)]
            mf = driver(mol)
            mf.xc = self.xc
            mf._numint.libxc = xcfun
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
        elif mf and (step == "tddft"):
            mf = pyscf.tddft.TDDFT(mf)
            mf.nstates = self.nstates
        elif mf and (step == "tda"):
            mf = pyscf.tddft.TDA(mf)
            mf.nstates = self.nstates
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
        mol.verbose = 4
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
        if self.root:
            grad_driver.state = self.root
        gradient = grad_driver.kernel()

        try:
            e_tot = mf._scf.e_tot
        except AttributeError:
            e_tot = mf.e_tot

        results = {
            "energy": e_tot,
            "forces": -gradient.flatten(),
        }

        self.mf = mf
        if self.track:
            if self.track_root(atoms, coords):
                # Redo the calculation with the updated root
                results = self.get_forces(atoms, coords)
        return results

    def get_hessian(self, atoms, coords):
        mol = self.prepare_input(atoms, coords)
        mf = self.run(mol)
        H = mf.Hessian().kernel()

        # The returned hessian is 4d ... ok. This probably serves a purpose
        # that I don't understand. We transform H to a nice, simple 2d array.
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
        self.calc_number += 1
        return mf

    def prepare_overlap_data(self):
        gs_mf = self.mf._scf
        exc_mf = self.mf

        if self.wfow == None:
            occ_num, virt_num = exc_mf.xy[0][0].shape
            self.wfow = WFOWrapper(occ_num, virt_num, calc_number=self.calc_number,
                                   basis=None, charge=None, out_dir=self.out_dir)
        gs_energy = gs_mf.e_tot
        mo_coeffs = gs_mf.mo_coeff

        fake_mos_fn = Path(self.make_fn("mos"))
        if not fake_mos_fn.exists():
            fake_mos_str = self.wfow.fake_turbo_mos(mo_coeffs)
            with open(fake_mos_fn, "w") as handle:
                handle.write(fake_mos_str)
        else:
            self.log("Skipping creation of MOs in TURBOMOLE format, as the file "
                     f"'{fake_mos_fn}' already exists.")

        # Shape = (nstates, 2 (X,Y), occ, virt)
        ci_coeffs = np.array(exc_mf.xy)
        # In TDA calculations Y is just the integer 0.
        if ci_coeffs.ndim == 2:
            ci_coeffs = np.array([state[0] for state in exc_mf.xy])
        # In TD-DFT calculations the Y vectors is also present
        else:
            # Sum X and Y to X+Y
            ci_coeffs = ci_coeffs.sum(axis=1)
        # Norm (X+Y) to 1 for every state
        norms = np.linalg.norm(ci_coeffs, axis=(1,2))
        ci_coeffs /= norms[:,None,None]

        exc_energies = exc_mf.e_tot
        all_energies = np.zeros(exc_energies.size + 1)
        all_energies[0] = gs_energy
        all_energies[1:] = exc_energies
        return fake_mos_fn, mo_coeffs, ci_coeffs, all_energies

    def __str__(self):
        return f"PySCF({self.name})"
