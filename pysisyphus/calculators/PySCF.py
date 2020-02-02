#!/usr/bin/env python3

import os
from pathlib import Path
import shutil

import numpy as np
import pyscf
from pyscf import gto, grad, lib, hessian, tddft, qmmm
from pyscf.dft import xcfun
# from pyscf.lib.chkfile import save_mol
# from pyscf.tools import cubegen

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
        "tda": ("dft", "tda"),
    }

    def __init__(self, basis, xc=None, method="scf",  mem=2000,
                 root=None, nstates=None, auxbasis=None, keep_chk=True,
                 **kwargs):
        super().__init__(**kwargs)

        self.basis = basis
        self.xc = xc
        self.method = method.lower()
        if self.xc and self.method != "tddft":
            self.method = "dft"
        # self.multistep = self.method in ("mp2 tddft".split())
        self.mem = mem
        self.root = root
        self.nstates = nstates
        if self.method == "tddft":
            assert self.nstates, "nstates must be set with method='tddft'!"
        if self.track:
            assert self.root <= self.nstates, "'root' must be smaller " \
                "than 'nstates'!"
        self.auxbasis = auxbasis
        self.keep_chk = keep_chk

        self.chkfile = None
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
            # TODO: grids
            # grids = pyscf.dft.gen_grid.Grids(mol)
            # grids.level = 4
            # grids.build()
            # mf.grids = grids
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
            raise Exception("Unknown method '{step}'!")
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
        # Personally, I patched mole.py so it doesn't print
        # messages regarding the output-file for verbose > QUIET.
        # Just uncomment the lines after
        #   if self.verbose > logger.QUIET:
        #       ...
        # in 'mole.Mole.build'. Around line 2046 for pyscf 1.6.5.
        mol.output = self.make_fn(self.out_fn)
        mol.max_memory = self.mem * self.pal
        mol.build(parse_arg=False)

        return mol

    def get_energy(self, atoms, coords, prepare_kwargs=None):
        if prepare_kwargs is None:
            prepare_kwargs = {}
        point_charges = prepare_kwargs.get("point_charges", None)

        mol = self.prepare_input(atoms, coords)
        mf = self.run(mol, point_charges=point_charges)
        results = {
            "energy": mf.e_tot,
        }

        return results

    def get_forces(self, atoms, coords, prepare_kwargs=None):
        if prepare_kwargs is None:
            prepare_kwargs = {}
        point_charges = prepare_kwargs.get("point_charges", None)

        mol = self.prepare_input(atoms, coords)
        mf = self.run(mol, point_charges=point_charges)
        # >>> mf.chkfile = '/path/to/chkfile'
        # >>> mf.init_guess = 'chkfile'
        grad_driver = mf.Gradients()
        if self.root:
            grad_driver.state = self.root
        gradient = grad_driver.kernel()
        self.log(f"Completed gradient step")

        try:
            e_tot = mf._scf.e_tot
        except AttributeError:
            e_tot = mf.e_tot

        results = {
            "energy": e_tot,
            "forces": -gradient.flatten(),
        }

        if self.track:
            self.store_overlap_data(atoms, coords)
            if self.track_root():
                # Redo the calculation with the updated root
                results = self.get_forces(atoms, coords)
        return results

    def get_hessian(self, atoms, coords, prepare_kwargs=None):
        if prepare_kwargs is None:
            prepare_kwargs = {}
        point_charges = prepare_kwargs.get("point_charges", None)

        mol = self.prepare_input(atoms, coords)
        mf = self.run(mol, point_charges=point_charges)
        H = mf.Hessian().kernel()

        # The returned hessian is 4d ... ok. This probably serves a purpose
        # that I don't understand. We transform H to a nice, simple 2d array.
        H = np.hstack(np.concatenate(H, axis=1))
        results = {
            "energy": mf.e_tot,
            "hessian": H,
        }

        return results

    def run(self, mol, point_charges=None):
        steps = self.multisteps[self.method]
        self.log(f"Running steps '{steps}' for method {self.method}")
        for i, step in enumerate(steps):
            if i == 0:
                mf = self.get_driver(step, mol=mol)
                assert step in ("scf", "dft")
                if self.chkfile:
                    # Copy old chkfile to new chkfile
                    new_chkfile = self.make_fn("chkfile", return_str=True)
                    shutil.copy(self.chkfile, new_chkfile)
                    self.chkfile = new_chkfile
                    mf.chkfile = self.chkfile
                    mf.init_guess = "chkfile"
                    self.log(f"Using '{self.chkfile}' as initial guess for {step} calculation.")
                if self.auxbasis:
                    mf.density_fit(auxbasis=self.auxbasis)
                    self.log(f"Using density fitting with auxbasis {self.auxbasis}.")

                if point_charges is not None:
                    mf = qmmm.mm_charge(mf, point_charges[:,:3], point_charges[:,3])
                    self.log(f"Added {len(point_charges)} point charges with "
                             f"sum(q)={sum(point_charges[:,3]):.4f}.")
            else:
                mf = self.get_driver(step, mf=prev_mf)  # noqa: F821

            if self.keep_chk and (self.chkfile is None) and (step in ("dft", "scf")):
                self.chkfile = self.make_fn("chkfile", return_str=True)
                try:
                    os.remove(self.chkfile)
                except FileNotFoundError:
                    self.log(f"Tried to remove '{self.chkfile}'. It doesn't exist.")
                self.log(f"Created chkfile '{self.chkfile}'")
                mf.chkfile = self.chkfile
            mf.kernel()
            self.log(f"Completed {step} step")
            prev_mf = mf

        # Keep mf and dump mol
        # save_mol(mol, self.make_fn("mol.chk"))
        self.mf = mf
        self.calc_counter += 1

        return mf

    def prepare_overlap_data(self):
        gs_mf = self.mf._scf
        exc_mf = self.mf

        gs_energy = gs_mf.e_tot
        mo_coeffs = gs_mf.mo_coeff.T

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
        return mo_coeffs, ci_coeffs, all_energies

    def parse_charges(self):
        # Mulliken charges
        results = self.mf.analyze(with_meta_lowdin=False)

        return results[0][1]

    def __str__(self):
        return f"PySCF({self.name})"
