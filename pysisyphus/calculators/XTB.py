#!/usr/bin/env python3

import glob
import os

import numpy as np
import pyparsing as pp

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.config import Config
from pysisyphus.constants import BOHR2ANG
from pysisyphus.calculators.parser import parse_turbo_gradient

from pysisyphus.xyzloader import make_xyz_str


class XTB(Calculator):

    def __init__(self, gbsa="", **kwargs):
        super(XTB, self).__init__(**kwargs)

        self.gbsa = gbsa
        self.uhf = self.mult - 1

        self.inp_fn = "xtb.xyz"
        self.out_fn = "xtb.out"
        self.to_keep = ("out", "gradient")

        self.parser_funcs = {
            "grad": self.parse_gradient,
        }

        self.base_cmd = Config["xtb"]["cmd"]

    def reattach(self, last_calc_cycle):
        pass

    def prepare_coords(self, atoms, coords):
        coords = coords * BOHR2ANG
        return make_xyz_str(atoms, coords.reshape((-1, 3)))

    def get_forces(self, atoms, coords):
        inp = self.prepare_coords(atoms, coords)
        add_args = f"-gfn -chrg {self.charge} -uhf {self.uhf} -grad".split()
        # Use solvent model if specified
        if self.gbsa:
            gbsa = f"-gbsa {self.gbsa}".split()
            add_args = add_args + gbsa
        results = self.run(inp, calc="grad", add_args=add_args)
        return results

    def parse_gradient(self, path):
        return parse_turbo_gradient(path)

    def __str__(self):
        return "XTB calculator"


if __name__ == "__main__":
    path = "/scratch/test/xtbtest"
    xtb = XTB()
    xtb.parse_gradient(path)
