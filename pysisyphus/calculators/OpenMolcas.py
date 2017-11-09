#!/usr/bin/env python3

import glob
import os
import re

import numpy as np
import pyparsing as pp

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.config import Config
from pysisyphus.constants import BOHR2ANG

from qchelper.geometry import make_xyz_str


class OpenMolcas(Calculator):

    def __init__(self):
        super(OpenMolcas, self).__init__()

        self.inp_fn = "openmolcas.in"
        self.out_fn = "openmolcas.out"

        self.openmolcas_input = """
        >> copy /scratch/test/ommin/excrp.es_opt.RasOrb  $Project.RasOrb
        &gateway
         coord
          {xyz_str}
         basis
          6-31G*
         group
          nosym
 
        &seward

        &rasscf
         charge
          1
         spin
          1
         fileorb
          $Project.RasOrb
         thrs
          1.0e-6,1.0e-2,1.0e-2
         ciroot
          2 2 1
         rlxroot
          2

        &alaska
         pnew
        """

        self.parser_funcs = {
            "grad": self.parse_gradient,
        }

        self.base_cmd = Config["openmolcas"]["cmd"]

    def prepare_coords(self, atoms, coords):
        coords = coords * BOHR2ANG
        return make_xyz_str(atoms, coords.reshape((-1, 3)))

    def get_forces(self, atoms, coords):
        xyz_str = self.prepare_coords(atoms, coords)
        inp = self.openmolcas_input.format(xyz_str=xyz_str)
        add_args = ("-clean", "-oe", self.out_fn)
        results = self.run(inp, calc="grad", add_args=add_args)
        return results

    def parse_gradient(self, path):
        results = {}
        gradient_fn = os.path.join(path, self.out_fn)
        with open(gradient_fn) as handle:
            text = handle.read()

        # Search for the block containing the gradient table
        regex = "Molecular gradients(.+?)--- Stop Module:\s*alaska"
        float_regex = "([\d\.\-]+)"
        floats = [float_regex for i in range(3)]
        line_regex = "([A-Z\d]+)\s*" + "\s*".join(floats)
        energy_regex = "RASSCF state energy =\s*" + float_regex
        energy = float(re.search(energy_regex, text).groups()[0])

        mobj = re.search(regex, text, re.DOTALL)
        gradient = list()
        for line in mobj.groups()[0].split("\n"):
            # Now look for the lines containing the gradient
            mobj = re.match(line_regex, line.strip())
            if not mobj:
                continue
            # Discard first column (atom+number)
            gradient.append(mobj.groups()[1:])
        gradient = np.array(gradient, dtype=np.float).flatten()

        results["energy"] = energy
        results["forces"] = -gradient

        return results

    def __str__(self):
        return "OpenMolcas calculator"


if __name__ == "__main__":
    path = "/scratch/test/openmolcas_calculator"
    om = OpenMolcas()
    #om.parse_gradient(path)
    from pysisyphus.helpers import geom_from_library
    geom = geom_from_library("dieniminium_cation_s1_opt.xyz")
    geom.set_calculator(om)
    forces = geom.forces
    print(forces)
