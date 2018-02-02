#!/usr/bin/env python3

import glob
import logging
import os
from pathlib import Path
import re
import shutil

import numpy as np
import pyparsing as pp

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import AU2EV
from pysisyphus.config import Config
from pysisyphus.calculators.parser import parse_turbo_gradient
from pysisyphus.calculators.WFOWrapper import WFOWrapper


class Turbomole(Calculator):

    def __init__(self, control_path, root=None,
                 track=False, **kwargs):
        super(Turbomole, self).__init__(**kwargs)

        self.control_path = Path(control_path)
        self.root = root
        self.track = track
        self.to_keep = ("control", "mos", "alpha", "beta", "out",
                        "ciss_a", "ucis_a")

        self.parser_funcs = {
            "force": self.parse_force,
        }

        self.inp_fn = ""
        self.out_fn = "turbomole.out"
        # MO coefficient files
        self.mos = None
        self.alpha = None
        self.beta = None

        # Prepare base_cmd
        with open(self.control_path / "control") as handle:
            text = handle.read()
        base_cmd = ["dscf", "grad"]
        # Check for RI
        if ("$rij" in text) or ("$rik" in text):
            base_cmd = ["ridft", "rdgrad"]
            self.log("Found RI calculation.")

        self.uhf = "$uhf" in text
        assert self.uhf == False, "Implement occ for UHF!"

        self.td = False
        # Check for excited state calculation
        if "$exopt" in text:
            exopt_re = "\$exopt\s*(\d+)"
            self.root = int(re.search(exopt_re, text)[1])
            base_cmd[1] = "egrad"
            self.prepare_td(text)
        if self.track:
            assert self.td, "track=True can't be used without $exopt!"
        # Right now this can't handle a root flip from some excited state
        # to the ground state ... Then we would need grad/rdgrad again,
        # instead of egrad.
        self.base_cmd = ";".join(base_cmd)
        self.log(f"Using base_cmd {self.base_cmd}")

    def prepare_td(self, text):
        self.log("Preparing for excited state gradient calculations")
        self.td = True
        # Determine number of occupied orbitals
        occ_re = "closed shells\s+(\w)\s*\d+-(\d+)"
        self.occ_mos = int(re.search(occ_re, text)[2])
        self.log(f"Found {self.occ_mos} occupied MOs.")
        # Number of spherical basis functions. May be different from CAO
        nbf_re = "nbf\(AO\)=(\d+)"
        nbf = int(re.search(nbf_re, text)[1])
        # Determine number of virtual orbitals
        self.virt_mos = nbf - self.occ_mos
        self.wfow = WFOWrapper(self.occ_mos, self.virt_mos)
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
        if calc_type != "force":
            raise Exception("Can only do force now.")
            """To rectify this we have to construct the basecmd
            dynamically and construct it ad hoc. We could set a RI flag
            in the beginning and select the correct scf binary here from
            it. Then we select the following binary on demand, e.g. aoforce
            or rdgrad or egrad etc."""
        path = self.prepare_path(use_in_run=True)
        # Copy everything from the reference control_dir into this path
        all_src_paths = self.control_path.glob("./*")
        """Maybe we shouldn't copy everything because it may give convergence
        problems?"""
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
        # Use inp=None because we got no special input...
        # Use shell=True because we have to chain commands like ridft;rdgrad
        results = self.run(None, calc="force", shell=True)
        if self.track:
            self.check_for_root_flip(atoms, coords)
        return results

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

    def ci_coeffs_above_thresh(self, eigenpair, thresh=1e-5):
        arr = np.array(eigenpair["vector"])
        arr = arr.reshape(self.occ_mos, -1)
        mo_inds = np.where(np.abs(arr) > thresh)
        #ci_coeffs = arr[mo_inds]
        #vals_sq = vals**2
        #for from_mo, to_mo, vsq, v in zip(*inds, vals_sq, vals):
        #    print(f"\t{from_mo+1} -> {to_mo+1+self.occ_mos} {v:.04f} {vsq:.02f}")
        #return ci_coeffs, mo_inds
        return arr, mo_inds

    def check_for_root_flip(self, atoms, coords):
        with open(self.td_vec_fn) as handle:
            text = handle.read()
        # Parse the eigenvectors
        eigenpair_list = self.parse_td_vectors(text)
        # Filter for relevant MO indices and corresponding CI coefficients
        coeffs_inds = [self.ci_coeffs_above_thresh(ep)
                       for ep in eigenpair_list]
        ci_coeffs, mo_inds = zip(*coeffs_inds)
        ci_coeffs = np.array(ci_coeffs)
        self.wfow.store_iteration(atoms, coords, self.mos, ci_coeffs, mo_inds)
        if self.calc_counter >= 1:
            new_root = self.wfow.track()
            if new_root != self.root:
                self.log("Found a root flip from {self.root} to {new_root}")
                self.root = new_root

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
            self.td_vec_fn = kept_fns["ciss_a"]


    def __str__(self):
        return "Turbomole calculator"
