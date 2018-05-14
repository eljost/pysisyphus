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

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import AU2EV
from pysisyphus.config import Config
from pysisyphus.calculators.parser import (parse_turbo_gradient,
                                           parse_turbo_ccre0_ascii,)
from pysisyphus.calculators.WFOWrapper import WFOWrapper


class Turbomole(Calculator):

    def __init__(self, control_path, root=None,
                 track=False, double_mol_path=None, **kwargs):
        super(Turbomole, self).__init__(**kwargs)

        self.control_path = Path(control_path)
        self.root = root
        self.track = track
        self.double_mol_path = Path(double_mol_path)

        # Check if the overlap matrix will be printed and assert
        # that no SCF iterations are done.
        if self.double_mol_path:
            with open(self.double_mol_path / "control") as handle:
                text = handle.read()
            assert (re.search("\$intsdebug\s*sao", text) and
                    re.search("\$scfiterlimit\s*0", text)), "Please set " \
                   "$intsdebug sao and $scfiterlimit 0 !"

        self.to_keep = ("control", "mos", "alpha", "beta", "out",
                        "ciss_a", "ucis_a",
                        "__ccre*", "exstates")

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
        self.wfow = WFOWrapper(self.occ_mos, self.virt_mos,
                               basis=None, charge=None,
                               calc_number=self.calc_number)
        self.td_vec_fn = None
        self.ci_coeffs = None
        self.mo_inds = None

    def prepare_coords(self, atoms, coords):
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
        # Write coordinates
        coord_str = self.prepare_coords(atoms, coords)
        coord_fn = path / "coord"
        with open(coord_fn, "w") as handle:
            handle.write(coord_str)
        # Copy MO coefficients from previous cycle with this calculator
        # if present.
        if self.mos:
            shutil.copy(self.mos, path / "mos")
            self.log(f"Using {self.mos}")
        elif self.alpha and self.beta:
            shutil.copy(self.alpha, path / "alpha")
            shutil.copy(self.beta, path / "beta")
            self.log(f"Using {self.alpha} and {self.beta}")

    def get_forces(self, atoms, coords):
        self.prepare_input(atoms, coords, "force")
        kwargs = {
                "calc": "force",
                "shell": True, # To allow 'ridft; rdgrad' etc.
                "hold": self.track, # Keep the files for WFOverlap
        }
        # Use inp=None because we don't use any dedicated input besides
        # the previously prepared control file and the current coords.
        results = self.run(None, **kwargs)
        if self.track:
            if self.track_root(atoms, coords):
                # Redo the calculation with the updated root
                results = self.get_forces(atoms, coords)
            self.calc_counter += 1
        return results

    def run_calculation(self, atoms, coords):
        """Basically some kind of dummy method that can be called
        to execute Turbomole with the stored cmd of this calculator."""
        self.prepare_input(atoms, coords, "noparse")
        kwargs = {
                "calc": "noparse",
                "shell": True,
                "hold": self.track, # Keep the files for WFOverlap
        }
        results = self.run(None, **kwargs)
        if self.track:
            self.track_root(atoms, coords)
            self.calc_counter += 1
        return results

    def run_double_mol_calculation(self, atoms, coords1, coords2):
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
        }
        results = self.run(None, **kwargs)
        self.calc_counter -= 1
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

    def parse_force(self, path):
        return parse_turbo_gradient(path)

    def parse_td_vectors(self, text):
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

    def store_wfo_data(self, atoms, coords):
        # Parse eigenvectors from escf/egrad calculation
        if self.second_cmd != "ricc2":
            with open(self.td_vec_fn) as handle:
                text = handle.read()
            eigenpair_list = self.parse_td_vectors(text)
            eigenpair_list = [ep["vector"] for ep in eigenpair_list]
        # Parse eigenvectors from ricc2 calculation
        else:
            eigenpair_list = [self.parse_cc2_vectors(ccre)
                              for ccre in self.ccres]
        self.wfow.store_iteration(atoms, coords, self.mos, eigenpair_list)

    def track_root(self, atoms, coords):
        """Store the information of the current iteration and if possible
        calculate the overlap with the previous iteration."""
        self.store_wfo_data(atoms, coords)
        # In the first iteration we have nothing to compare to
        old_root = self.root
        if self.calc_counter >= 1:
            last_two_coords = self.wfow.last_two_coords
            ao_ovlp = self.run_double_mol_calculation(atoms, *last_two_coords)
            self.root = self.wfow.track(old_root=self.root, ao_ovlp=ao_ovlp)
            if self.root != old_root:
                self.log("Found a root flip from {old_root} to {self.root}!")

        # True if a root flip occured
        return not (self.root == old_root)

    def keep(self, path):
        kept_fns = super().keep(path)
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
                self.td_vec_fn = kept_fns["ciss_a"]
            elif self.ricc2:
                self.ccres = kept_fns["ccres"]
            else:
                raise Exception("Something went wrong!")

    def run_after(self, path):
        # Convert binary CCRE0 files to ASCII for easier parsing
        for ccre in path.glob("CCRE0-*"):
            cmd = f"ricctools -dump {ccre.name}".split()
            result = subprocess.Popen(cmd, cwd=path,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
            result.wait()

    def __str__(self):
        return "Turbomole calculator"
