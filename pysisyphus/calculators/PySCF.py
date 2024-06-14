import os
import shutil
from typing import Optional

import numpy as np
import pyscf
from pyscf import gto, grad, lib, hessian, tddft, qmmm

from pysisyphus.calculators.OverlapCalculator import OverlapCalculator
from pysisyphus.elem_data import ATOMIC_NUMBERS
from pysisyphus.helpers import geom_loader
from pysisyphus.wavefunction import Shell, Wavefunction
from pysisyphus.wavefunction.shells_pyscf import PySCFShells
from pysisyphus.wavefunction.helpers import BFType


def from_pyscf_mol(mol, **kwargs):
    shells = list()
    for bas_id in range(mol.nbas):
        L = mol.bas_angular(bas_id)
        center = mol.bas_coord(bas_id)
        coeffs = mol.bas_ctr_coeff(bas_id).flatten()
        exps = mol.bas_exp(bas_id)
        assert coeffs.size == exps.size, "General contractions are not yet supported."
        center_ind = mol.bas_atom(bas_id)
        atom_symbol = mol.atom_symbol(center_ind)
        atomic_num = ATOMIC_NUMBERS[atom_symbol.lower()]
        shell = Shell(L, center, coeffs, exps, center_ind, atomic_num)
        shells.append(shell)
    return PySCFShells(shells, **kwargs)


def kernel_to_wavefunction(mf):
    mol = mf.mol
    shells = from_pyscf_mol(mol)
    wf = Wavefunction(
        atoms=tuple(shells.atoms),
        coords=shells.coords3d,
        charge=mol.charge,
        mult=mol.multiplicity,
        unrestricted=mf.mo_occ.ndim == 2,
        occ=mol.nelec,
        C=mf.mo_coeff,
        bf_type=BFType.CARTESIAN if mol.cart else BFType.PURE_SPHERICAL,
        shells=shells,
        mo_ens=mf.mo_energy,
    )
    return wf


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
        "scf": ("scf",),
        "tdhf": ("scf", "tddft"),
        "tdahf": ("scf", "tda"),
        "dft": ("dft",),
        "mp2": ("scf", "mp2"),
        "tddft": ("dft", "tddft"),
        "tda": ("dft", "tda"),
    }
    pruning_method = {
        "nwchem": pyscf.dft.gen_grid.nwchem_prune,
        "sg1": pyscf.dft.gen_grid.sg1_prune,
        "treutler": pyscf.dft.gen_grid.treutler_prune,
        "none": None,
    }

    def __init__(
        self,
        basis,
        xc=None,
        method="scf",
        auxbasis=None,
        keep_chk=True,
        verbose: int = 0,
        unrestricted: Optional[bool] = None,
        grid_level: int = 3,
        pruning="nwchem",
        td_triplets=False,
        **kwargs,
    ):
        super().__init__(**kwargs)

        self.basis = basis
        self.xc = xc
        self.method = method.lower()
        self.do_es = self.method in ("tda", "tddft", "tdhf", "tdahf")
        # Set method to DFT when an XC-functional was selected, but DFT wasn't explicitly
        # enabled. When an ES method was chosen it is kept.
        if self.xc and self.method not in ("tddft", "tda"):
            self.method = "dft"
        if self.do_es:
            assert self.nroots, f"nroots must be set with method='{self.method}'!"
        self.auxbasis = auxbasis
        self.keep_chk = keep_chk
        self.verbose = int(verbose)
        if unrestricted is None:
            unrestricted = self.mult > 1
        self.unrestricted = unrestricted
        self.grid_level = int(grid_level)
        self.pruning = pruning.lower()
        self.td_triplets = td_triplets

        self.chkfile = None
        self.out_fn = "pyscf.out"

        lib.num_threads(self.pal)

    @staticmethod
    def geom_from_fn(fn, **kwargs):
        geom = geom_loader(fn)
        geom.set_calculator(PySCF(**kwargs))
        return geom

    def set_scf_params(self, mf):
        mf.conv_tol = 1e-8
        mf.max_cycle = 150

    def build_grid(self, mf):
        mf.grids.level = self.grid_level
        mf.grids.prune = self.pruning_method[self.pruning]
        mf.grids.build()

    def prepare_mf(self, mf):
        # Method can be overriden in a subclass to modify the mf object.
        return mf

    def configure_td_mf(self, mf):
        mf.nstates = self.nroots
        # Whether to calculate spin-adapated triplets from a singlet reference
        # When triplets are requested, mf.singlet must be set to False.
        mf.singlet = not self.td_triplets

    def get_driver(self, step, mol=None, mf=None):
        def _get_driver():
            return self.drivers[(step, self.unrestricted)]

        if mol and (step == "dft"):
            driver = _get_driver()
            mf = driver(mol)
            mf.xc = self.xc
            self.set_scf_params(mf)
            self.build_grid(mf)
            mf = self.prepare_mf(mf)
        elif mol and (step == "scf"):
            driver = _get_driver()
            mf = driver(mol)
            self.set_scf_params(mf)
            mf = self.prepare_mf(mf)
        elif mf and (step == "mp2"):
            mp2_mf = _get_driver()
            mf = mp2_mf(mf)
        elif mf and (step == "tddft"):
            mf = pyscf.tddft.TDDFT(mf)
            self.configure_td_mf(mf)
        elif mf and (step == "tda"):
            mf = pyscf.tddft.TDA(mf)
            self.configure_td_mf(mf)
        else:
            raise Exception("Unknown method '{step}'!")
        return mf

    def prepare_mol(self, atoms, coords, build=True):
        mol = gto.Mole()
        mol.atom = [(atom, c) for atom, c in zip(atoms, coords.reshape(-1, 3))]
        mol.basis = self.basis
        mol.unit = "Bohr"
        mol.charge = self.charge
        mol.spin = self.mult - 1
        mol.symmetry = False
        mol.verbose = self.verbose
        # Personally, I patched mole.py so it doesn't print
        # messages regarding the output-file for verbose > QUIET.
        # Just uncomment the lines after
        #   if self.verbose > logger.QUIET:
        #       ...
        # in 'mole.Mole.build'. Around line 2046 for pyscf 1.6.5.
        # Search for "output file" in gto/mole.py
        # Search for "Large deviations found" in scf/{uhf,dhf,ghf}.py
        mol.output = self.make_fn(self.out_fn)
        mol.max_memory = self.mem * self.pal
        if build:
            mol.build(parse_arg=False)
        return mol

    def prepare_input(self, atoms, coords, build=True):
        mol = self.prepare_mol(atoms, coords, build=build)
        assert mol._built, "Please call mol.build(parse_arg=False)!"
        return mol

    def store_and_track(self, results, func, atoms, coords, **prepare_kwargs):
        if self.track:
            self.store_overlap_data(atoms, coords)
            if self.track_root():
                # Redo the calculation with the updated root
                results = func(atoms, coords, **prepare_kwargs)
        results["all_energies"] = self.parse_all_energies()
        return results

    def get_energy(self, atoms, coords, **prepare_kwargs):
        point_charges = prepare_kwargs.get("point_charges", None)

        mol = self.prepare_input(atoms, coords)
        mf = self.run(mol, point_charges=point_charges)
        all_energies = self.parse_all_energies()
        if self.root is not None:
            energy = all_energies[self.root]
        else:
            energy = all_energies[0]

        results = {
            "energy": energy,
        }
        results = self.store_and_track(
            results, self.get_energy, atoms, coords, **prepare_kwargs
        )
        return results

    def get_all_energies(self, atoms, coords, **prepare_kwargs):
        results = self.get_energy(atoms, coords, **prepare_kwargs)
        all_energies = self.parse_all_energies()
        results["all_energies"] = all_energies
        return results

    def get_forces(self, atoms, coords, **prepare_kwargs):
        point_charges = prepare_kwargs.get("point_charges", None)

        mol = self.prepare_input(atoms, coords)
        mf = self.run(mol, point_charges=point_charges)
        grad_driver = mf.Gradients()
        if self.root is not None:
            grad_driver.state = self.root
        gradient = grad_driver.kernel()
        self.log("Completed gradient step")

        if self.root is not None:
            # Index 0 of e_tot contains energy of 1st excited state etc..., so we have
            # to substract 1.
            e_tot = mf.e_tot[self.root - 1]
        else:
            e_tot = mf._scf.e_tot

        results = {
            "energy": e_tot,
            "forces": -gradient.flatten(),
        }
        results = self.store_and_track(
            results, self.get_forces, atoms, coords, **prepare_kwargs
        )
        return results

    def get_hessian(self, atoms, coords, **prepare_kwargs):
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
        # results = self.store_and_track(
        # results, self.get_hessian, atoms, coords, **prepare_kwargs
        # )
        return results

    def run_calculation(self, atoms, coords, **prepare_kwargs):
        return self.get_energy(atoms, coords, **prepare_kwargs)

    def get_wavefunction(self, atoms, coords, **prepare_kwargs):
        point_charges = prepare_kwargs.get("point_charges", None)

        method = self.multisteps[self.method][0]
        mol = self.prepare_input(atoms, coords)
        mf = self.run(mol, point_charges=point_charges, method=method)
        all_energies = self.parse_all_energies()
        energy = all_energies[0]
        results = {
            "energy": energy,
            "wavefunction": kernel_to_wavefunction(mf),
        }
        return results

    def run(self, mol, point_charges=None, method=None):
        if method is None:
            method = self.method
        steps = self.multisteps[method]
        self.log(f"Running steps '{steps}' for method {method}")
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
                    self.log(
                        f"Using '{self.chkfile}' as initial guess for {step} calculation."
                    )
                if self.auxbasis:
                    mf.density_fit(auxbasis=self.auxbasis)
                    self.log(f"Using density fitting with auxbasis {self.auxbasis}.")

                if point_charges is not None:
                    mf = qmmm.mm_charge(mf, point_charges[:, :3], point_charges[:, 3])
                    self.log(
                        f"Added {len(point_charges)} point charges with "
                        f"sum(q)={sum(point_charges[:,3]):.4f}."
                    )
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

    def parse_all_energies(self, exc_mf=None):
        if exc_mf is None:
            exc_mf = self.mf

        try:
            gs_energy = exc_mf._scf.e_tot
            exc_energies = exc_mf.e_tot
            all_energies = np.zeros(exc_energies.size + 1)
            all_energies[0] = gs_energy
            all_energies[1:] = exc_energies
        except AttributeError:
            gs_energy = exc_mf.e_tot
            all_energies = np.array((gs_energy,))

        return all_energies

    def prepare_overlap_data(self, path):
        gs_mf = self.mf._scf
        exc_mf = self.mf

        C = gs_mf.mo_coeff

        first_Y = exc_mf.xy[0][1]
        # In TDA calculations Y is just the integer 0.
        if isinstance(first_Y, int) and (first_Y == 0):
            X = np.array([state[0] for state in exc_mf.xy])
            Y = np.zeros_like(X)
        # In TD-DFT calculations the Y vectors is also present
        else:
            # Shape = (nstates, 2 (X,Y), occ, virt)
            ci_coeffs = np.array(exc_mf.xy)
            X = ci_coeffs[:, 0]
            Y = ci_coeffs[:, 1]

        all_energies = self.parse_all_energies(exc_mf)
        return C, X, Y, all_energies

    def prepare_overlap_data(self, path):
        gs_mf = self.mf._scf
        exc_mf = self.mf

        # MO coefficients
        if self.unrestricted:
            Ca, Cb = gs_mf.mo_coeff
        else:
            Ca = gs_mf.mo_coeff
            Cb = Ca.copy()

        xy = exc_mf.xy
        X = np.array([x for x, _ in xy])
        Y = np.array([y for _, y in xy])

        if self.unrestricted:
            Xa = X[:, 0]
            Xb = X[:, 1]
            Ya = Y[:, 0]
            Yb = Y[:, 1]
        else:
            Xa = X
            Xb = Xa.copy()
            Ya = Y
            Yb = Ya.copy()

        # Methods that don't yield a Y-vector (TDA/TDAHF) will only produce a 0
        # instread of a proper Y vector. Below we fix this by creating appropriate
        # arrays.
        if self.method in ("tda", "tdahf"):
            Ya = np.zeros_like(Xa)
            Yb = np.zeros_like(Xb)

        all_energies = self.parse_all_energies(exc_mf)
        return Ca, Xa, Ya, Cb, Xb, Yb, all_energies

    def parse_charges(self):
        results = self.mf.analyze(with_meta_lowdin=False)
        # Mulliken charges
        charges = results[0][1]
        return charges

    def get_chkfiles(self):
        return {
            "chkfile": self.chkfile,
        }

    def set_chkfiles(self, chkfiles):
        try:
            chkfile = chkfiles["chkfile"]
            self.chkfile = chkfile
            self.log(f"Set chkfile '{chkfile}' as chkfile.")
        except KeyError:
            self.log("Found no chkfile information in chkfiles!")

    def __str__(self):
        return f"PySCF({self.name})"
