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
O          0.00000        0.00000        0.11779
H          0.00000        0.75545       -0.47116
H          0.00000       -0.75545       -0.47116
*
"""

class ORCA(Calculator):

    def __init__(self): 
        super(ORCA, self).__init__()

        self.inp_fn = "orca.inp"
        self.out_fn = "orca.out"
        self.base_cmd = "/scratch/programme/orca_4_0_0_2_linux_x86-64/orca"

    def get_energy(self, coords):
        raise Exception("Not implemented!")

    def get_forces(self, coords):
        results = self.run(inp)
        return results

    def get_hessian(self, coords):
        raise Exception("Not implemented!")

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
            force = np.array(engrad[:3*atoms], dtype=np.float)
            results["energy"] = energy
            results["force"] = force

        return results

    def __str__(self):
        return "ORCA calculator"

if __name__ == "__main__":
    orca = ORCA()
    orca.parse("/tmp/tmpoqcz2lce")
