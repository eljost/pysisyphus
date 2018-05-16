#!/usr/bin/env python3

from collections import namedtuple
import logging
import os
from pathlib import Path
import re
import subprocess

import numpy as np
import pyparsing as pp

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.calculators.WFOWrapper import WFOWrapper
from pysisyphus.constants import AU2EV
from pysisyphus.config import Config


class Gaussian16(Calculator):

    def __init__(self, method, basis, nstates=None, root=None, **kwargs):
        super(Gaussian16, self).__init__(**kwargs)

        self.method = method
        self.basis = basis
        self.nstates = nstates
        self.root = root
        # When root or nstates is set, the other option is required too!
        if root or nstates:
            assert (root and nstates), "nstates and root have to "\
                                       "be given together!"
            self.wfow = None

        self.to_keep = ("fchk", "log", "dump_635r")

        self.fn_base = "gaussian16"
        self.inp_fn = f"{self.fn_base}.com"
        self.out_fn = f"{self.fn_base}.log"
        self.chk_fn = f"{self.fn_base}.chk"
        self.dump_base_fn = f"{self.fn_base}_rwfdump"

        self.gaussian_input = """
        %chk={chk_fn}
        #P {calc_type} {method}/{basis} {td} iop(9/40=4)
        # nosymm density=all pop=full
        
        title

        {charge} {mult}
        {coords}



        """

        self.parser_funcs = {
            "force": self.parse_force,
        }

        self.base_cmd = Config["gaussian16"]["cmd"]
        self.formchk_cmd = Config["gaussian16"]["formchk_cmd"]

    def make_td_str(self):
        if not self.root:
            return ""
        td = f"td=(nstates={self.nstates},root={self.root})"
        return td

    def prepare_input(self, atoms, coords, calc_type):
        coords = self.prepare_coords(atoms, coords)
        inp = self.gaussian_input.format(
                        chk_fn=self.chk_fn,
                        calc_type=calc_type,
                        method=self.method,
                        basis=self.basis,
                        charge=self.charge,
                        mult=self.mult,
                        coords=coords,
                        td=self.make_td_str(),
        )
        return inp

    def make_fchk(self, path):
        cmd = f"{self.formchk_cmd} {self.chk_fn}".split()
        result = subprocess.run(cmd, stdout=subprocess.PIPE, cwd=path)
        self.log("Created .fchk")

    def run_rwfdump(self, path, rwf_index):
        chk_path = path / self.chk_fn
        dump_fn = path / f"{self.dump_base_fn}_dump_{rwf_index}"
        cmd = f"rwfdump {chk_path} {dump_fn} {rwf_index}".split()
        proc = subprocess.run(cmd)
        self.log(f"Dumped {rwf_index} from .chk.")
        return dump_fn

    def run_after(self, path):
        # Create the .fchk file so we can keep it and parse it later on.
        self.make_fchk(path)
        if self.root:
            self.run_rwfdump(path, "635r")
            self.nmos, self.roots = self.parse_log(path)

    def parse_fchk(self, fchk_path, keys):
        with open(fchk_path) as handle:
            text = handle.read()

        def to_float(s, loc, toks):
            return float(toks[0])

        # Matches -4.10693837E-16 and 1.60184209E-15
        float_ = pp.Word(pp.nums + "+-E.").setParseAction(to_float)
        # Start with Empty so we can progessively build the
        # parser for all keys.
        parser = pp.Empty()
        def parser_for_key(key):
            return pp.Group(pp.Suppress(pp.SkipTo(key)) + key + pp.restOfLine
                            + pp.ZeroOrMore(float_))
        for key in keys:
            parser += parser_for_key(key)
        results = parser.parseString(text)
        results_dict = {}
        for key, res in zip(keys, results):
            # This handles scalar entries like
            # Total Energy  [...] R     -9.219940072302333E+01
            if len(res) == 2:
                results_dict[key] = float(res[-1].split()[-1])
            # This handles matrices like
            # Cartesian Gradient [...] R   N=           9 \
            # [Matrix entries]
            if len(res) > 2:
                results_dict[key] = np.array(res[2:])
        return results_dict

    def get_forces(self, atoms, coords):
        inp = self.prepare_input(atoms, coords, "force")
        kwargs = {
            "calc": "force",
            "hold": True,
        }
        results = self.run(inp, **kwargs)
        if self.root:
            self.store_wfo_data(atoms, coords)
        return results

    def parse_tddft(self, path):
        with open(path / self.out_fn) as handle:
            text = handle.read()
        td_re = "Excited State\s*\d+:\s*[\w\?-]+\s*([\d\.-]+?)\s*eV"
        matches = re.findall(td_re, text)
        assert len(matches) == self.nstates
        # Excitation energies in eV
        exc_energies = np.array(matches, dtype=np.float)
        # Convert to Hartree
        exc_energies /= AU2EV
        return exc_energies

    def parse_log(self, path):
        self.log(f"Parsing {self.out_fn}")

        def parse(text, regex, func):
            mobj = re.search(regex, text)
            return func(mobj[1])

        log_path = path / self.out_fn
        with open(log_path) as handle:
            text = handle.read()

        roots_re = "Root\s+(\d+)"
        roots = np.array(re.findall(roots_re, text), dtype=int).max()

        # NBasis=    16 NAE=    12 NBE=    12 NFC=     6 NFV=     0
        basis_re = "NBasis=(.+)NAE=(.+)NBE=(.+)NFC=(.+)NFV=(.+)"
        basis_mobj = re.search(basis_re, text)
        basis_funcs, alpha, beta, _, _ = [int(n) for n in basis_mobj.groups()]
        a_occ = alpha
        b_occ = beta
        a_vir = basis_funcs - a_occ
        b_vir = basis_funcs - b_occ
        restricted = alpha == beta
        act_re = "NROrb=(.*)NOA=(.*)NOB=(.*)NVA=(.*)NVB=(.*)"
        act_mobj = re.search(act_re, text)
        _, a_act, b_act, _, _ = [int(n) for n in act_mobj.groups()]
        a_core = a_occ - a_act
        b_core = b_occ - b_act

        NMOs = namedtuple("NMOs", "a_core, a_occ a_act a_vir "
                                  "b_core b_occ b_act b_vir "
                                  "restricted")

        nmos = NMOs(a_core, a_occ, a_act, a_vir,
                    b_core, b_occ, b_act, b_vir,
                    restricted)
        self.log(str(nmos))
        return nmos, roots

    def parse_635r_dump(self, dump_path, roots, nmos):
        self.log(f"Parsing 635r dump '{dump_path}'")
        with open(dump_path) as handle:
            text = handle.read()
        regex = "read left to right\):\s*(.+)"
        mobj = re.search(regex, text, re.DOTALL)
        arr_str = mobj[1].replace("D", "E")
        # Drop the first 12 items as they are always 0
        tmp = np.array(arr_str.split()[12:], dtype=np.float64)

        # the core electrons are frozen in TDDFT/TDA
        expected = (nmos.a_act*nmos.a_vir + nmos.b_act*nmos.b_vir) * roots * 2
        print(f"Expecting {expected} items. There are {tmp.size} elements.")
        coeffs = tmp[:expected]
        # 1. dim: X+Y, X-Y -> 2
        # 2. dim: roots -> variable
        # 3. dim: alpha, beta -> 2
        # the rest depends on the number of alpha and beta electrons
        if nmos.restricted:
            XpY, XmY = coeffs.reshape(2, roots, 2, -1)
            X = 0.5 * (XpY + XmY)
        # In the unrestricted case we can't use 2 for the 3. dimension anymore
        # (same number of alpha and beta electrons). The sizes for the respective
        # spins would be (a_act*a_vir) and (b_act*b_vir).
        else:
            raise Exception("Unrestricted not supported yet!")

        # xpy_fn = f"{dump_fn}_XpY"
        # np.save(xpy_fn, XpY)
        # xmy_fn = f"{dump_fn}_XmY"
        # np.save(xmy_fn, XmY)
        # x_fn = f"{dump_fn}_X"
        # np.save(x_fn, X)

        # Drop duplicate entries
        if self.nmos.restricted:
            # Drop beta part and restrict to the requested states
            X = X[:self.nstates, 0, :]

        X_full = np.zeros((self.nstates, nmos.a_occ, nmos.a_vir))
        X_full[:, nmos.a_act:] = X.reshape(-1, nmos.a_act, nmos.a_vir)

        return X_full

    def store_wfo_data(self, atoms, coords):
        # Create the WFOWrapper object if it is not already there
        if self.wfow == None:
            assert (self.nmos.restricted)
            occ_num, virt_num = self.nmos.a_occ, self.nmos.a_vir
            self.wfow = WFOWrapper(occ_num, virt_num, calc_number=self.calc_number,
                                   basis=None, charge=None)
        # Parse X eigenvector from 635r dump
        eigenpair_list = self.parse_635r_dump(self.dump_635r, self.roots, self.nmos)
        print(eigenpair_list)
        # Parse mo coefficients from .fchk file and write a 'fake' turbomole
        # mos file.
        keys = ("Alpha Orbital Energies", "Alpha MO coefficients")
        energies_key, mo_key = keys
        fchk_dict = self.parse_fchk(self.fchk, keys=keys)
        mo_coeffs = fchk_dict[mo_key]
        mo_energies = fchk_dict[energies_key]
        mo_coeffs = mo_coeffs.reshape(-1, mo_energies.size)
        fake_mos_str = self.wfow.fake_turbo_mos(mo_coeffs)
        fake_mos_fn = self.out_dir / self.make_fn("mos")
        with open(fake_mos_fn, "w") as handle:
            handle.write(fake_mos_str)
        self.wfow.store_iteration(atoms, coords, fake_mos_fn, eigenpair_list)

    def parse_force(self, path):
        results = {}
        keys = ("Total Energy", "Cartesian Gradient")
        fchk_path = Path(path) / f"{self.fn_base}.fchk"
        fchk_dict = self.parse_fchk(fchk_path, keys)
        results["energy"] = fchk_dict["Total Energy"]
        results["forces"] = -fchk_dict["Cartesian Gradient"]

        if self.nstates:
            # This sets the proper excited state energy in the
            # results dict and also stores all energies.
            exc_energies = self.parse_tddft(path)
            # G16 root input is 1 based, so we substract 1 to get
            # the right index here.
            root_exc_en = exc_energies[self.root-1]
            gs_energy = results["energy"]
            # Add excitation energy to ground state energy.
            results["energy"] += root_exc_en
            # Create a new array including the ground state energy
            # to save the energies of all states.
            all_ens = np.full(len(exc_energies)+1, gs_energy)
            all_ens[1:] += exc_energies
            results["tddft_energies"] = all_ens

        return results

        """
        # Parse the .log
        path = Path(path) / self.out_fn
        with open(path) as handle:
            text = handle.read()
        force_re = "\(Quad\)   \(Total\)(.+?)Item"
        mobj = re.search(force_re, text, re.DOTALL)
        force_table = mobj[1].strip()
        force = [line.strip().split()[2] for line in force_table.split("\n")]
        force = np.array(force, dtype=np.float)
        energy_re = "SCF Done:\s*E\(.+?\)\s*=\s*([\d\.\-]+)\s*A\.U."
        mobj = re.search(energy_re, text)
        energy = float(mobj[1])
        results["energy"] = energy
        results["forces"] = force
        return results
        """

    def keep(self, path):
        kept_fns = super().keep(path)
        self.fchk = kept_fns["fchk"]
        if self.root:
            self.dump_635r = kept_fns["dump_635r"]

    def __str__(self):
        return "Gaussian16 calculator"


if __name__ == "__main__":
    from pysisyphus.helpers import geom_from_library
    """
    geom = geom_from_library("dieniminium_cation_s1_opt.xyz")
    charge = 1
    mult = 1
    root = 1
    nstates = 1
    """
    geom = geom_from_library("hcn.xyz")
    charge = 0
    mult = 1
    method = "b3lyp"
    basis = "sto-3g"
    g16 = Gaussian16(method, basis, charge=charge, mult=mult)
    geom.set_calculator(g16)
    #forces = geom.forces
    #print(forces)
    # Total SCF Density
    # Total CI Rho(1) Density
    g16.nstates = 2
    g16.root = 1
    path = Path("/scratch/programme/pysisyphus/tests/test_g16_but2en_iso/crashed_image_0")
    g16.parse_force(path)

