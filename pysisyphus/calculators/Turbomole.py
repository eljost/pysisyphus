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


class Turbomole(Calculator):

    def __init__(self, control_path, root=None, **kwargs):
        super(Turbomole, self).__init__(**kwargs)

        self.control_path = Path(control_path)
        self.root = root
        self.to_keep = ("mos", "alpha", "beta", "out")

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
        if "$exopt" in text:
            base_cmd[1] = "egrad"
            self.log("Using excited state gradients.")
        self.uhf = "$uhf" in text
        # Right now this can't handle a root flip from some excited state
        # to the ground state ... Then we would need grad/rdgrad again,
        # instead of egrad.
        self.base_cmd = ";".join(base_cmd)
        self.log(f"Using base_cmd {self.base_cmd}")

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
        return results

    def parse_force(self, path):
        return parse_turbo_gradient(path)

    def keep(self, path):
        kept_fns = super().keep(path)
        if self.uhf:
            self.alpha = kept_fns["alpha"]
            self.beta = kept_fns["beta"]
        else:
            self.mos = kept_fns["mos"]
        # Maybe copy more files like the vectors from egrad
        # sing_a, trip_a, dipl_a etc.

    def __str__(self):
        return "Turbomole calculator"


if __name__ == "__main__":
    from pysisyphus.helpers import geom_from_library
    geom = geom_from_library("h2o.xyz")
    control_path = Path("/scratch/wfoverlap_1.0/pyscf/h2o_pure")
    turbo = Turbomole(control_path)
    atoms, coords = geom.atoms, geom.coords
    #turbo.prepare_input(atoms, coords, "inp_type")
    #coord_str = turbo.prepare_coords(atoms, coords)
    #print(coord_str)
    #geom.set_calculator(turbo)
    #forces = geom.forces
    #print(forces)
    #p = Path("/tmp/calculator_0_000_11vo_teg")
    #turbo.parse_force(p)

    fn = "/tmp/calculator_0_000_ab7q7o9y"
    path = Path(fn)
    turbo.keep(path)
