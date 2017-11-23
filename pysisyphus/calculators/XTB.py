#!/usr/bin/env python3

import glob
import os

import numpy as np
import pyparsing as pp

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.config import Config
from pysisyphus.constants import BOHR2ANG

from pysisyphus.xyzloader import make_xyz_str


class XTB(Calculator):

    def __init__(self, **kwargs):
        super(XTB, self).__init__(**kwargs)

        self.uhf = self.mult - 1

        self.inp_fn = "xtb.xyz"
        self.out_fn = "xtb.out"

        self.parser_funcs = {
            "grad": self.parse_gradient,
        }

        self.base_cmd = Config["xtb"]["cmd"]

    def prepare_coords(self, atoms, coords):
        coords = coords * BOHR2ANG
        return make_xyz_str(atoms, coords.reshape((-1, 3)))

    def get_forces(self, atoms, coords):
        inp = self.prepare_coords(atoms, coords)
        add_args = f"-gfn -chrg {self.charge} -uhf {self.uhf} -grad".split()
        results = self.run(inp, calc="grad", add_args=add_args)
        return results

    def parse_gradient(self, path):
        results = {}
        gradient_fn = glob.glob(os.path.join(path, "gradient"))
        if not gradient_fn:
            raise Exception("XTB gradient file not found!")
        assert(len(gradient_fn) == 1)
        gradient_fn = gradient_fn[0]
        with open(gradient_fn) as handle:
            text = handle.read()

        def to_float(s, loc, toks):
            match = toks[0].replace("D", "E")
            return float(match)
        float_ = pp.Word(pp.nums + ".-D+").setParseAction(to_float)

        cycle = pp.Word(pp.nums).setResultsName("cycle")
        scf_energy = float_.setResultsName("scf_energy")
        grad_norm = float_.setResultsName("grad_norm")
        float_line = float_ + float_ + float_
        coord_line = pp.Group(float_line + pp.Word(pp.alphas))
        grad_line = pp.Group(float_line)

        parser = (
            pp.Literal("$grad") +
            pp.Literal("cycle =") + cycle +
            pp.Literal("SCF energy =") + scf_energy +
            pp.Literal("|dE/dxyz| =") + grad_norm +
            pp.OneOrMore(coord_line).setResultsName("coords") +
            pp.OneOrMore(grad_line).setResultsName("grad") +
            pp.Literal("$end")
        )
        parsed = parser.parseString(text)
        gradient = np.array(parsed["grad"].asList()).flatten()

        results["energy"] = parsed["scf_energy"]
        results["forces"] = -gradient

        return results

    def __str__(self):
        return "XTB calculator"


if __name__ == "__main__":
    path = "/scratch/test/xtbtest"
    xtb = XTB()
    xtb.parse_gradient(path)
