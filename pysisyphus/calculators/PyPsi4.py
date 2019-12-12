#!/usr/bin/env python3

try:
    import psi4
except ImportError:
    pass
import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import BOHR2ANG
from pysisyphus.xyzloader import make_xyz_str

class PyPsi4(Calculator):

    def __init__(self, method, basis, mem=2000, **kwargs):
        super().__init__(**kwargs)

        self.method = method
        self.basis = basis
        self.mem = mem

        self.meth_bas = f"{self.method}/{self.basis}"

        psi4.core.be_quiet()
        psi4.set_memory(f"{mem} MB")

    def get_forces(self, atoms, coords):
        xyz = make_xyz_str(atoms, coords.reshape(-1,3)*BOHR2ANG)
        mol = psi4.geometry(xyz)
        energy, wfn = psi4.energy(self.meth_bas, return_wfn=True)
        gradient = psi4.gradient(self.meth_bas, molecule=mol, ref_wfn=wfn)

        results = {
            "energy": energy,
            "forces": -np.asarray(gradient).flatten(),
        }
        return results

    def __str__(self):
        return f"PyPsi4({self.name})"
