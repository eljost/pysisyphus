#!/usr/bin/env python3

import glob
import logging
import os
from pathlib import Path
import re
import shutil
import struct

import numpy as np
import pyparsing as pp

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import AU2EV
from pysisyphus.config import Config
from pysisyphus.calculators.parser import parse_turbo_gradient
from pysisyphus.calculators.WFOWrapper import WFOWrapper


class Turbomole(Calculator):

    def __init__(self, control_path, root=None,
                 track=False, wfo_basis=None, **kwargs):
        super(Turbomole, self).__init__(**kwargs)

        self.control_path = Path(control_path)
        self.root = root
        self.track = track
        self.wfo_basis = wfo_basis
        self.to_keep = ("control", "mos", "alpha", "beta", "out",
                        "ciss_a", "ucis_a", "__CCRE*")

        self.parser_funcs = {
            "force": self.parse_force,
            "tddft": self.parse_tddft,
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
        if self.track:
            assert (self.td or self.ricc2), "track=True can only be used " \
                "in connection with excited state calculations."
            assert self.wfo_basis != None, "Please supply a valid basis " \
                "supported by pyscf, to calculate the AO overlaps."
        # Right now this can't handle a root flip from some excited state
        # to the ground state ... Then we would need grad/rdgrad again,
        # instead of egrad.
        self.scf_cmd = scf_cmd
        self.second_cmd = second_cmd
        self.base_cmd = ";".join((self.scf_cmd, self.second_cmd))
        self.log(f"Using base_cmd {self.base_cmd}")

    def prepare_td(self, text):
        self.log("Preparing for excited state (gradient) calculations")
        self.wfow = WFOWrapper(self.occ_mos, self.virt_mos, self.wfo_basis,
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
        if calc_type not in ("force", "tddft"):
            raise Exception("Can only do force and tddft for now.")
            """To rectify this we have to construct the basecmd
            dynamically and construct it ad hoc. We could set a RI flag
            in the beginning and select the correct scf binary here from
            it. Then we select the following binary on demand, e.g. aoforce
            or rdgrad or egrad etc."""
        path = self.prepare_path(use_in_run=True)
        # Copy everything from the reference control_dir into this path
        all_src_paths = self.control_path.glob("./*")
        """Maybe we shouldn't copy everything because it may give convergence
        problems? Right now we use the initial MO guess generated in the
        reference path for all images along the path."""
        globs = [p for p in all_src_paths]
        for glob in self.control_path.glob("./*"):
            shutil.copy(glob, path)
        # Write coordinates
        coord_str = self.prepare_coords(atoms, coords)
        coord_fn = path / "coord"
        with open(coord_fn, "w") as handle:
            handle.write(coord_str)
        # Copy MO coefficients from previous calculation if present
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
            if self.check_for_root_flip(atoms, coords):
                # Redo the calculation with the updated root
                results = self.get_forces(atoms, coords)
            self.calc_counter += 1
        return results

    def get_tddft(self, atoms, coords):
        """See get_forces"""
        self.prepare_input(atoms, coords, "tddft")
        kwargs = {
                "calc": "tddft",
                "shell": True,
                "hold": self.track, # Keep the files for WFOverlap
        }
        results = self.run(None, **kwargs)
        if self.track:
            self.check_for_root_flip(atoms, coords)
            self.calc_counter += 1
        return results

    def parse_force(self, path):
        return parse_turbo_gradient(path)

    def parse_tddft(self, path):
        return None

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
        """Adapted from TheoDORE 1.6 lib/file_parser.py
        https://sourceforge.net/projects/theodore-qc/
        
        Maybe there is a mistake in TheoDORE as only 20 bytes are skipped
        before reading the coeffs. Maybe we have to skip 28 bytes?!"""

        with open(ccre, "rb") as handle:
            handle.read(8)
            method = struct.unpack('8s', handle.read(8))[0]
            handle.read(8)
            nentry = struct.unpack('i', handle.read(4))[0]
            handle.read(28) 

            assert(nentry % self.virt_mos == 0)
            active_mos = nentry // self.virt_mos
            frozen_mos = self.occ_mos - active_mos

            coeffs = [struct.unpack("d", handle.read(8))[0]
                      for _ in range(active_mos*self.virt_mos)]
            eigenpairs = np.array(coeffs).reshape(-1, self.virt_mos)
            # 'eigenpairs' has the shape (active occ. MOs x virtual MOs)
            # and doesn't include the frozen orbitals.
            # In the end we need an array with the shape
            #  (occ. MOs x virt. MOs), that is including frozen orbitals.
            # Right now we are basically missing some columns, so we create
            # the full array and add the previous 'eigenpairs' block.
            eigenpairs_full = np.zeros((self.occ_mos, self.virt_mos))
            eigenpairs_full[frozen_mos:,:] = eigenpairs
        return eigenpairs_full

    def ci_coeffs_above_thresh(self, eigenpair, thresh=1e-5):
        arr = np.array(eigenpair)
        arr = arr.reshape(self.occ_mos, -1)
        mo_inds = np.where(np.abs(arr) > thresh)
        #ci_coeffs = arr[mo_inds]
        #vals_sq = vals**2
        #for from_mo, to_mo, vsq, v in zip(*inds, vals_sq, vals):
        #    print(f"\t{from_mo+1} -> {to_mo+1+self.occ_mos} {v:.04f} {vsq:.02f}")
        #return ci_coeffs, mo_inds
        return arr, mo_inds

    def check_for_root_flip(self, atoms, coords):
        """Call WFOverlap, store the information of the current iteration and
        calculate the overlap with the previous iteration, if possible."""

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
        # Filter for relevant MO indices and corresponding CI coefficients
        coeffs_inds = [self.ci_coeffs_above_thresh(ep)
                       for ep in eigenpair_list]
        ci_coeffs, mo_inds = zip(*coeffs_inds)
        ci_coeffs = np.array(ci_coeffs)
        self.wfow.store_iteration(atoms, coords, self.mos, ci_coeffs, mo_inds)
        # In the first iteration we have nothing to compare to
        old_root = self.root
        if self.calc_counter >= 1:
            self.root = self.wfow.track(old_root=self.root)
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

    def __str__(self):
        return "Turbomole calculator"
