#!/usr/bin/env python3

try:
    from ase.atoms import Atoms
except ImportError:
    pass

import sys
try:
    xtb_path = "/scratch/programme/xtb-190806/python"
    sys.path.append(xtb_path)
    from xtb_mod import GFN2
except ImportError:
    pass


from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import BOHR2ANG


class PyXTB(Calculator):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.gfn2 = GFN2()

    def prepare_input(self, atoms, coords, calc_type):
        return None

    def calculate(self, atoms, coords):
        atoms = Atoms(atoms, coords.reshape(-1, 3)*BOHR2ANG)
        self.gfn2.calculate(atoms)
        res = self.gfn2.results
        return res

    def get_forces(self, atoms, coords):
        res = self.calculate(atoms, coords)
        results = {
            "forces": res["forces"].flatten(),
            "energy": res["energy"],
        }
        return results

    def __str__(self):
        return "PyXTB calculator"
