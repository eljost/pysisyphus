#!/usr/bin/env python3

import numpy as np
import qcengine as qcng
import qcelemental as qcel

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import BOHR2ANG


class QCEngine(Calculator):

    def __init__(self, program, model, keywords=None, **kwargs):
        super().__init__(**kwargs)

        self.program = program
        self.model = model
        self.keywords = keywords

        # TODO: pal, memory

    def get_molecule(self, atoms, coords):
        c3d = coords.reshape(-1, 3)# * BOHR2ANG
        mol_inp = {
            "symbols": atoms,
            "geometry": c3d,
            "molecular_charge": self.charge,
            "molecular_multiplicity": self.mult,
        }
        molecule = qcel.models.Molecule.from_data(mol_inp)

        return molecule

    def compute(self, molecule, driver):
        inp = {
            "molecule": molecule,
            "driver": driver,
            "model": self.model,
            "keywords": self.keywords,
        }

        result = qcng.compute(inp, program=self.program, return_dict=True)
        # with open("stdout", "w") as handle:
            # handle.write(__["stdout"])
        return result

    def get_energy(self, atoms, coords, prepare_kwargs=None):
        mol = self.get_molecule(atoms, coords)

        res = self.compute(mol, driver="energy")

        results = {
            "energy": res["return_result"],
        }

        return results

    def get_forces(self, atoms, coords, prepare_kwargs=None):
        mol = self.get_molecule(atoms, coords)

        res = self.compute(mol, driver="gradient")

        results = {
            "energy": res["properties"]["return_energy"],
            "forces": -np.array(res["return_result"]),
        }

        return results

    def __str__(self):
        return f"QCECalculator({self.name})"
