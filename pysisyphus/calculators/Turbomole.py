#!/usr/bin/env python3

import glob
import logging
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess

import numpy as np
import pyparsing as pp

from pysisyphus.calculators.OverlapCalculator import OverlapCalculator
from pysisyphus.constants import AU2EV
from pysisyphus.calculators.parser import (parse_turbo_gradient,
                                           parse_turbo_ccre0_ascii,
                                           parse_turbo_mos,
                                           parse_turbo_exstates)


class Turbomole(OverlapCalculator):

    conf_key = "turbomole"

    def __init__(self, control_path, root=None,
                 double_mol_path=None, mem=1000, **kwargs):
        super(Turbomole, self).__init__(**kwargs)

        self.control_path = Path(control_path)
        self.root = root
        self.double_mol_path = double_mol_path
        if self.double_mol_path:
            self.double_mol_path = Path(self.double_mol_path)
        self.mem = int(mem)

        # Check if the overlap matrix will be printed and assert
        # that no SCF iterations are done.
        if self.double_mol_path:
            with open(self.double_mol_path / "control") as handle:
                text = handle.read()
            assert (re.search("\$intsdebug\s*sao", text) and
                    re.search("\$scfiterlimit\s*0", text)), "Please set " \
                   "$intsdebug sao and $scfiterlimit 0 !"

        self.to_keep = ("control", "mos", "alpha", "beta", "out",
                        "ciss_a", "ucis_a", "gradient", "sing_a",
                        "__ccre*", "exstates", "coord",
                        "mwfn_wf:wavefunction.molden",
                        "input.xyz",
        )

        self.parser_funcs = {
            "force": self.parse_force,
            "double_mol": self.parse_double_mol,
            "noparse": lambda path: None,
        }

        # Turbomole uses the 'control' file implicitly
        self.inp_fn = ""
        self.out_fn = "turbomole.out"
        # MO coefficient files
        self.mos = None
        self.alpha = None
        self.beta = None

        # Prepare base_cmd
        with open(self.control_path / "control") as handle:
            text = handle.read()
        scf_cmd = "dscf"
        second_cmd = "grad"
        # Check for RI
        if ("$rij" in text) or ("$rik" in text):
            scf_cmd = "ridft"
            second_cmd = "rdgrad"
            self.log("Found RI calculation.")

        self.uhf = "$uhf" in text
        assert self.uhf == False, "Implement occ for UHF!"

        # Determine number of occupied orbitals
        occ_re = "closed shells\s+(\w)\s*\d+-(\d+)"
        self.occ_mos = int(re.search(occ_re, text)[2])
        self.log(f"Found {self.occ_mos} occupied MOs.")
        # Number of spherical basis functions. May be different from CAO
        nbf_re = "nbf\(AO\)=(\d+)"
        nbf = int(re.search(nbf_re, text)[1])
        # Determine number of virtual orbitals
        self.virt_mos = nbf - self.occ_mos

        self.td = False
        self.ricc2 = False
        # Check for excited state calculation
        if "$exopt" in text:
            exopt_re = "\$exopt\s*(\d+)"
            self.root = int(re.search(exopt_re, text)[1])
            second_cmd = "egrad"
            self.prepare_td(text)
            self.td = True
        elif "$soes" in text:
            second_cmd = "escf"
            self.td = True
            self.prepare_td(text)
        elif ("$ricc2" in text) and ("$excitations" in text):
            self.ricc2 = True
            second_cmd = "ricc2"
            self.prepare_td(text)
            self.root = self.get_ricc2_root(text)
            self.frozen_mos = int(re.search("implicit core=\s*(\d+)", text)[1])
            self.log(f"Found {self.frozen_mos} frozen orbitals.")
        if self.track:
            assert (self.td or self.ricc2), "track=True can only be used " \
                "in connection with excited state calculations."
        # Right now this can't handle a root flip from some excited state
        # to the ground state ... Then we would need grad/rdgrad again,
        # instead of egrad.
        self.scf_cmd = scf_cmd
        self.second_cmd = second_cmd
        self.base_cmd = ";".join((self.scf_cmd, self.second_cmd))
        self.log(f"Using base_cmd {self.base_cmd}")

    def get_ricc2_root(self, text):
        regex = "geoopt.+?state=\((.+?)\)"
        mobj = re.search(regex, text)
        if not mobj:
            root = None
        elif mobj[1] == "x":
            root = 0
        else:
            assert mobj[1].startswith("a "), "symmetry isn't supported!"
            root = int(mobj[1][-1])
        return root

    def prepare_td(self, text):
        self.log("Preparing for excited state (gradient) calculations")
        self.td_vec_fn = None
        self.ci_coeffs = None
        self.mo_inds = None

    def prepare_turbo_coords(self, atoms, coords):
        fmt = "{:<20.014f}"
        coord_str = "$coord\n"
        for atom, coord in zip(atoms, coords.reshape(-1, 3)):
            coord_line = (fmt+fmt+fmt).format(*coord) + atom.lower() + "\n"
            coord_str += coord_line
        coord_str += "$end"
        return coord_str

    def prepare_input(self, atoms, coords, calc_type):
        if calc_type not in ("force", "double_mol", "noparse"):
            raise Exception("Can only do force and tddft for now.")
            """To rectify this we have to construct the basecmd
            dynamically and construct it ad hoc. We could set a RI flag
            in the beginning and select the correct scf binary here from
            it. Then we select the following binary on demand, e.g. aoforce
            or rdgrad or egrad etc."""
        path = self.prepare_path(use_in_run=True)
        if calc_type == "double_mol":
            copy_from = self.double_mol_path
        else:
            copy_from = self.control_path
        # Copy everything from the reference control_dir into this path
        # Use self.control_path for all calculations except the double
        # molecule calculation.
        all_src_paths = copy_from.glob("./*")
        """Maybe we shouldn't copy everything because it may give convergence
        problems? Right now we use the initial MO guess generated in the
        reference path for all images along the path."""
        globs = [p for p in all_src_paths]
        for glob in copy_from.glob("./*"):
            shutil.copy(glob, path)
        xyz_str = self.prepare_xyz_string(atoms, coords)
        with open(path / "input.xyz", "w") as handle:
            handle.write(xyz_str)
        # Write coordinates
        coord_str = self.prepare_turbo_coords(atoms, coords)
        coord_fn = path / "coord"
        with open(coord_fn, "w") as handle:
            handle.write(coord_str)
        # Copy MO coefficients from previous cycle with this calculator
        # if present.
        if self.mos:
            shutil.copy(self.mos, path / "mos")
            self.log(f"Using {self.mos} as MO guess.")
        elif self.alpha and self.beta:
            shutil.copy(self.alpha, path / "alpha")
            shutil.copy(self.beta, path / "beta")
            self.log(f"Using {self.alpha} and {self.beta} as MO guesses.")
        if self.td_vec_fn:
            # The suffix contains the true name with a leading
            # dot, that we drop.
            td_vec_fn = self.td_vec_fn.suffix[1:]
            shutil.copy(self.td_vec_fn, path / td_vec_fn)
            self.log(f"Using '{self.td_vec_fn}' as escf guess.")
        root_log_msg = "with current root information."
        if self.root and self.td:
            repl = f"$exopt {self.root}"
            self.sub_control("\$exopt\s*(\d+)", f"$exopt {self.root}",
                             root_log_msg)
            self.log(f"Using '{repl}'")
        if self.root and self.ricc2:
            repl = f"state=(a {self.root})"
            log_msg = " with current root."
            self.sub_control("state=\(a\s+(?P<state>\d+)\)", f"state=(a {self.root})",
                             root_log_msg)
            self.log(f"Using '{repl}' for geoopt.")

    def sub_control(self, pattern, repl, log_msg="", **kwargs):
        path = self.path_already_prepared
        assert path
        self.log(f"Updating control file in '{path}' {log_msg}")
        control_path = path / "control"
        with open(control_path) as handle:
            text = handle.read()
        text = re.sub(pattern, repl, text, **kwargs)
        with open(control_path, "w") as handle:
            handle.write(text)

    def get_pal_env(self):
        env_copy = os.environ.copy()
        env_copy["PARA_ARCH"] = "SMP"
        env_copy["PARNODES"] = str(self.pal)
        env_copy["SMPCPUS"] = str(self.pal)

        return env_copy

    def get_energy(self, atoms, coords):
        results = self.get_forces(atoms, coords)
        del results["forces"]
        return results

    def get_forces(self, atoms, coords, cmd=None):
        self.prepare_input(atoms, coords, "force")
        kwargs = {
                "calc": "force",
                "shell": True, # To allow chained commands like 'ridft; rdgrad'
                "hold": self.track, # Keep the files for WFOverlap
                "env": self.get_pal_env(),
                "cmd": cmd,
        }
        # Use inp=None because we don't use any dedicated input besides
        # the previously prepared control file and the current coords.
        results = self.run(None, **kwargs)
        if self.track:
            prev_run_path = self.last_run_path
            self.store_overlap_data(atoms, coords)
            root_flipped = self.track_root()
            self.calc_counter += 1
            if root_flipped:
                # Redo gradient calculation for new root.
                results = self.get_forces(atoms, coords, cmd=self.second_cmd)
            self.last_run_path = prev_run_path
        shutil.rmtree(self.last_run_path)
        return results

    def run_calculation(self, atoms, coords):
        """Basically some kind of dummy method that can be called
        to execute Turbomole with the stored cmd of this calculator."""
        self.prepare_input(atoms, coords, "noparse")
        kwargs = {
                "calc": "noparse",
                "shell": True,
                # "hold": self.track, # Keep the files for WFOverlap
                "env": self.get_pal_env(),
        }
        results = self.run(None, **kwargs)
        if self.track:
            self.store_overlap_data(atoms, coords)
        return results

    def run_double_mol_calculation(self, atoms, coords1, coords2):
        if not self.double_mol_path:
            self.log("Skipping double molecule calculations as double mol "
                     "path is not specified.!")
            return None
        self.log("Running double molecule calculation")
        double_atoms = atoms + atoms
        double_coords = np.hstack((coords1, coords2))
        self.prepare_input(double_atoms, double_coords, "double_mol")
        kwargs = {
                "calc": "double_mol",
                "shell": True,
                "keep": False,
                "hold": True,
                "cmd": self.scf_cmd,
                "env": self.get_pal_env(),
        }
        results = self.run(None, **kwargs)
        return results

    def parse_double_mol(self, path):
        """Parse a double molecule overlap matrix from Turbomole output
        to be used with WFOWrapper."""
        with open(path / self.out_fn) as handle:
            text = handle.read()
        regex = "OVERLAP\(SAO\)\s+-+([\d\.E\-\s*\+]+)\s+-+"
        ovlp_str = re.search(regex, text)[1]
        ovlp = np.array(ovlp_str.strip().split(), dtype=np.float64)
        mo_num = self.occ_mos + self.virt_mos
        double_mo_num = 2 * mo_num
        full_ovlp = np.zeros((double_mo_num, double_mo_num))
        full_ovlp[np.tril_indices(double_mo_num)] = ovlp
        double_mol_S = full_ovlp[mo_num:,:mo_num]
        return double_mol_S

    def parse_mos(self):
        pass

    def parse_force(self, path):
        results = parse_turbo_gradient(path)
        return results

    def parse_td_vectors(self, text):
        """For TDA calculations only the X vector is present in the
        ciss_a/etc. file. In TDDFT calculations there are twise as
        much items compared with TDA. The first half corresponds to
        (X+Y) and the second half to (X-Y). X can be calculated as
        X = ((X+Y)+(X-Y))/2. Y is then given as Y = (X+Y)-X. The
        normalization can then by checked as
            np.concatenate((X, Y)).dot(np.concatenate((X, -Y)))
        and should be 1."""
        def to_float(s, loc, toks):
            match = toks[0].replace("D", "E")
            return float(match)

        float_ = pp.Word(pp.nums + ".-D+").setParseAction(to_float)
        integer = pp.Word(pp.nums).setParseAction(lambda t: int(t[0]))
        float_chrs = pp.nums + "D.+-"
        float_20 = pp.Word(float_chrs, exact=20).setParseAction(to_float)

        title = pp.Literal("$title")
        symmetry = pp.Literal("$symmetry") \
                   + pp.Word(pp.alphanums).setResultsName("symmetry")
        tensor_dim = pp.Literal("$tensor space dimension") \
                    + integer.setResultsName("tensor_dim")
        scfinstab = pp.Literal("$scfinstab") \
                    + pp.Word(pp.alphanums).setResultsName("scfinstab")
        subspace_dim = pp.Literal("$current subspace dimension") \
                       + integer.setResultsName("subspace_dim")
        converged = pp.Literal("$current iteration converged")
        eigenpairs = pp.Literal("$eigenpairs")
        eigenpair = pp.Group(
                        integer.setResultsName("state")
                        + pp.Literal("eigenvalue =")
                        + float_.setResultsName("eigenvalue")
                        + pp.Group(pp.OneOrMore(float_20)).setResultsName("vector"))
        end = pp.Literal("$end")

        parser = (title
                  + symmetry
                  + tensor_dim
                  + scfinstab
                  + subspace_dim
                  + converged
                  + eigenpairs
                  + pp.OneOrMore(eigenpair).setResultsName("eigenpairs")
                  + end
        )
        result = parser.parseString(text)
        states = result["subspace_dim"]
        eigenpairs = result["eigenpairs"]
        eigenpair_list = [eigenpairs[i].asDict() for i in range(states)]
        return eigenpair_list

    def parse_cc2_vectors(self, ccre):
        with open(ccre) as handle:
            text = handle.read()
        coeffs = parse_turbo_ccre0_ascii(text)
        coeffs = coeffs.reshape(-1, self.virt_mos)

        eigenpairs_full = np.zeros((self.occ_mos, self.virt_mos))
        eigenpairs_full[self.frozen_mos:,:] = coeffs
        from_inds, to_inds = np.where(np.abs(eigenpairs_full) > .1)

        # for i, (from_, to_) in enumerate(zip(from_inds, to_inds)):
            # sq = eigenpairs_full[from_, to_]**2
            # print(f"{from_+1:02d} -> {to_+self.occ_mos+1:02d}: {sq:.2%}")
        # print()

        return eigenpairs_full

    def parse_gs_energy(self):
        """Several places are possible:
            $subenergy from control file
            total energy from turbomole.out
            Final MP2 energy from turbomole.out with ADC(2)
            Final CC2 energy from turbomole.out with CC(2)
        """
        float_re = "([\d\-\.E]+)"
        regexs = [
                  # CC2 ground state energy
                  ("out", "Final CC2 energy\s*:\s*" + float_re, 0),
                  # ADC(2) ground state energy
                  ("out", "Final MP2 energy\s*:\s*" + float_re, 0),
                  ("control", "\$subenergy.*$\s*" + float_re, re.MULTILINE),
                  # DSCF ground state energy
                  ("out", "total energy\s*=\s*" + float_re, 0),
        ]
        for file_attr, regex, flag in regexs:
            regex_ = re.compile(regex, flags=flag)
            with open(getattr(self, file_attr)) as handle:
                text = handle.read()
            mobj = regex_.search(text)
            try:
                gs_energy = float(mobj[1])
                self.log(f"Parsed ground state energy from '{file_attr}' using "
                         f"regex '{regex[:11]}'.")
                return gs_energy
            except TypeError:
                continue
        raise Exception("Couldn't parse ground state energy!")

    def prepare_overlap_data(self):
        # Parse eigenvectors from escf/egrad calculation
        gs_energy = self.parse_gs_energy()
        if self.second_cmd != "ricc2":
            self.log(f"Reading CI coefficients from '{self.td_vec_fn}'.")
            with open(self.td_vec_fn) as handle:
                text = handle.read()
            ci_coeffs = self.parse_td_vectors(text)
            exc_energies = [cc["eigenvalue"] for cc in ci_coeffs]
            ci_coeffs = [cc["vector"] for cc in ci_coeffs]
            all_energies = np.full(len(exc_energies)+1, gs_energy)
            all_energies[1:] += exc_energies
            # s1 = np.array(ci_coeffs[0])
            # # XmY, XpY = s1.reshape(2, -1)
            # XpY, XmY = s1.reshape(2, -1)
            # X = (XmY+XpY)/2
            # Y = XpY-X
        # Parse eigenvectors from ricc2 calculation
        else:
            ci_coeffs = [self.parse_cc2_vectors(ccre)
                         for ccre in self.ccres]
            with open(self.exstates) as handle:
                exstates_text = handle.read()
            exc_energies_by_model = parse_turbo_exstates(exstates_text)
            # Drop CCS and take energies from whatever model was used
            exc_energies = [(model, exc_ens) for model, exc_ens in exc_energies_by_model
                            if model != "CCS"]
            assert len(exc_energies) == 1
            model, exc_energies = exc_energies[0]
            all_energies = np.full(len(exc_energies)+1, gs_energy)
            all_energies[1:] += exc_energies
            self.log("Parsing of all energies for ricc2 is not yet implemented!")
        ci_coeffs = np.array(ci_coeffs)
        states = ci_coeffs.shape[0]
        X_len = self.occ_mos * self.virt_mos
        if ci_coeffs.shape[1] == (2*X_len):
            self.log("TDDFT calculation with X and Y vectors present. "
            )
            ci_coeffs = ci_coeffs[:,:X_len]
        ci_coeffs = ci_coeffs.reshape(states, self.occ_mos, self.virt_mos)
        self.log(f"Reading MO coefficients from '{self.mos}'.")
        with open(self.mos) as handle:
            text = handle.read()
        mo_coeffs = parse_turbo_mos(text)
        self.log(f"Reading electronic energies from '{self.out}'.")
        return mo_coeffs, ci_coeffs, all_energies

    def keep(self, path):
        kept_fns = super().keep(path)
        self.out = kept_fns["out"]
        self.control = kept_fns["control"]
        if self.uhf:
            self.alpha = kept_fns["alpha"]
            self.beta = kept_fns["beta"]
        else:
            self.mos = kept_fns["mos"]
        # Maybe copy more files like the vectors from egrad
        # sing_a, trip_a, dipl_a etc.
        assert"ucis_a" not in kept_fns, "Implement for UKS TDA"
        if self.track:
            if self.td:
                td_key_present = [k for k in ("ciss_a", "sing_a", "ucis_a")
                                  if k in kept_fns][0]
                self.td_vec_fn = kept_fns[td_key_present]
            elif self.ricc2:
                self.ccres = kept_fns["ccres"]
                self.exstates = kept_fns["exstates"]
            else:
                raise Exception("Something went wrong!")
            self.mwfn_wf = kept_fns["mwfn_wf"]

    def run_after(self, path):
        # Convert binary CCRE0 files to ASCII for easier parsing
        for ccre in path.glob("CCRE0-*"):
            cmd = f"ricctools -dump {ccre.name}".split()
            result = subprocess.Popen(cmd, cwd=path,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
            result.wait()

        if self.td:
            self.make_molden(path)
        # With ricc2 we probably have a frozen core that we have to disable
        # temporarily before creating the molden file. Afterwards we restore
        # the original control file with the frozen core.
        elif self.ricc2:
            # Backup original control file
            ctrl_backup = path / "control.backup"
            shutil.copy(path / "control", ctrl_backup)
            # We have to remove line with implicit core in the control file
            with open(path / "control") as handle:
                text = handle.read()
            lines = text.split("\n")
            lines = [l for l in lines if "implicit core" not in l]
            with open(path / "control", "w") as handle:
                handle.write("\n".join(lines))
            self.make_molden(path)
            # Restore control backup
            shutil.copy(ctrl_backup, path / "control")

    def make_molden(self, path):
        cmd = "tm2molden norm".split()
        fn = "wavefunction.molden"
        stdin = f"""{fn}

        """
        res = subprocess.Popen(cmd, cwd=path, universal_newlines=True,
                               stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
        stdout, stderr = res.communicate(stdin)
        res.terminate()

    def propagate_wavefunction(self, calc):
        if self.mos:
            calc.mos = self.mos
        elif self.uhf and self.alpha and self.beta:
            calc.alpha = self.alpha
            calc.beta = self.beta

    def __str__(self):
        return "Turbomole calculator"
