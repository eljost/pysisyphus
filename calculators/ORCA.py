#!/usr/bin/env python3

import glob
import os
import re
import subprocess

import numpy as np

from calculators.Calculator import Calculator

inp="""!BP86 def2-SV(P) def2/J TightSCF engrad

%maxcore 1000

*xyz 0 1
{}
*
"""

class ORCA(Calculator):

    def __init__(self): 
        super(ORCA, self).__init__()

        self.inp_fn = "orca.inp"
        self.out_fn = "orca.out"
        self.base_cmd = "/scratch/programme/orca_4_0_0_2_linux_x86-64/orca"

    def get_energy(self, atoms, coords):
        raise Exception("Not implemented!")

    def get_forces(self, atoms, coords):
        coords = coords.reshape(-1, 3)
        coords = "\n".join(
            ["{} {} {} {}".format(a, *c) for a, c in zip(atoms, coords)]
        )
        #import logging; logging.info("Running calculation!")
        results = self.run(inp.format(coords))
        return results

    def parse(self, path):
        results = {}
        engrad_fn = glob.glob(os.path.join(path, "*.engrad"))
        if engrad_fn:
            assert(len(engrad_fn) == 1)
            engrad_fn = engrad_fn[0]
            with open(engrad_fn) as handle:
                engrad = handle.read()
            engrad = re.findall("([\d\-\.]+)", engrad)
            atoms = int(engrad.pop(0))
            energy = float(engrad.pop(0))
            force = -np.array(engrad[:3*atoms], dtype=np.float) * 0.529177249
            results["energy"] = energy
            results["forces"] = force

        return results

    def __str__(self):
        return "ORCA calculator"
