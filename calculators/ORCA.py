#!/usr/bin/env python3

import glob
import os
import re
import subprocess

import numpy as np
import pyparsing as pp

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.config import Config

#inp="""!BP86 def2-SV(P) def2/J TightSCF engrad
inp="""!HF def2-SV(P) TightSCF NoRI {}

*xyz 0 1
{}
*
"""

def make_sym_mat(table_block):
    mat_size = int(table_block[1])
    # Orca prints blocks of 5 columns
    arr = np.array(table_block[2:], dtype=np.float)
    assert(arr.size == mat_size**2)
    block_size = 5*mat_size
    cbs = [arr[i*block_size:(i+1)*block_size].reshape(mat_size, -1)
           for i in range(arr.size // block_size + 1)
    ]
    return np.concatenate(cbs, axis=1)

class ORCA(Calculator):

    def __init__(self): 
        super(ORCA, self).__init__()

        self.inp_fn = "orca.inp"
        self.out_fn = "orca.out"

        self.parser_funcs = {
            "grad": self.parse_engrad,
            "hessian": self.parse_hessian,
        }

        self.base_cmd = Config["orca"]["cmd"]

    def get_energy(self, atoms, coords):
        raise Exception("Not implemented!")

    def get_forces(self, atoms, coords):
        coords = coords.reshape(-1, 3)
        coords = "\n".join(
            ["{} {} {} {}".format(a, *c) for a, c in zip(atoms, coords)]
        )
        results = self.run(inp.format("engrad", coords), calc="grad")
        return results

    def get_hessian(self, atoms, coords):
        coords = coords.reshape(-1, 3)
        coords = "\n".join(
            ["{} {} {} {}".format(a, *c) for a, c in zip(atoms, coords)]
        )
        results = self.run(inp.format("freq", coords), calc="hessian")
        return results

    def parse_hessian(self, path):
        results = {}
        hessian_fn = glob.glob(os.path.join(path, "*.hess"))
        assert(len(hessian_fn) == 1)
        hessian_fn = hessian_fn[0]
        if not hessian_fn:
            raise Exception("ORCA calculation failed.")
        with open(hessian_fn) as handle:
            text = handle.read()

        integer = pp.Word(pp.nums)
        float_ = pp.Word(pp.nums + ".-")
        plus = pp.Literal("+")
        minus = pp.Literal("-")
        E = pp.Literal("E")
        scientific = pp.Combine(float_ + E + pp.Or([plus, minus]) + integer)

        table_header_line = pp.Suppress(integer + pp.restOfLine)
        scientific_line = pp.Suppress(integer) + pp.OneOrMore(scientific)
        scientific_block = table_header_line + pp.OneOrMore(scientific_line)
        float_line = pp.Suppress(integer) + float_
        comment_line = pp.Literal("#") + pp.restOfLine
        mass_xyz_line = (pp.Word(pp.alphas) + float_ +
                         pp.Group(pp.OneOrMore(float_))
        )

        block_name = pp.Word(pp.alphas + "$_")
        block_length = integer

        block_int = block_name + block_length
        block_float = block_name + float_
        block_table = block_name + integer + pp.OneOrMore(scientific_block)
        block_table_two_int = (block_name + integer + pp.Suppress(integer)
                               + pp.OneOrMore(scientific_block)
        )
        block_float_table = block_name + integer + pp.OneOrMore(float_line)
        block_atoms = block_name + integer + pp.OneOrMore(mass_xyz_line)

        act_atom = block_int.setResultsName("act_atom")
        act_coord = block_int.setResultsName("act_coord")
        act_energy = block_float.setResultsName("act_energy")
        hessian = block_table.setResultsName("hessian")
        vib_freqs = block_float_table.setResultsName("vib_freqs")
        normal_modes = block_table_two_int.setResultsName("normal_modes")
        atoms = block_atoms.setResultsName("atoms")

        parser = (block_name + act_atom + act_coord + act_energy
                  + hessian + vib_freqs + normal_modes
                  + pp.OneOrMore(comment_line) + atoms
        )
        parsed = parser.parseString(text)
        results["hessian"] = make_sym_mat(parsed["hessian"])
        return results

    def parse_engrad(self, path):
        results = {}
        engrad_fn = glob.glob(os.path.join(path, "*.engrad"))
        if not engrad_fn:
            raise Exception("ORCA calculation failed.")

        assert(len(engrad_fn) == 1)
        engrad_fn = engrad_fn[0]
        with open(engrad_fn) as handle:
            engrad = handle.read()
        engrad = re.findall("([\d\-\.]+)", engrad)
        atoms = int(engrad.pop(0))
        energy = float(engrad.pop(0))
        force = -np.array(engrad[:3*atoms], dtype=np.float) * 0.529177249
        import logging; logging.warning("still using anstroem internally for orca...")
        results["energy"] = energy
        results["forces"] = force

        return results

    def __str__(self):
        return "ORCA calculator"
