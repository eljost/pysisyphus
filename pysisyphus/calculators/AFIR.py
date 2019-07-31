#!/usr/bin/env python3

import itertools as it

import autograd
import autograd.numpy as anp
import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.elem_data import COVALENT_RADII


def afir_closure(fragment_indices, cov_radii, gamma, rho=1, p=6):
    """rho=1 pushes fragments together, rho=-1 pulls fragments apart."""

    # See https://onlinelibrary.wiley.com/doi/full/10.1002/qua.24757
    # Eq. (9) for extension to 3 fragments.
    assert len(fragment_indices) == 2

    inds = np.array(list(it.product(*fragment_indices)))
    cov_rad_sums = cov_radii[inds].sum(axis=1)

    # 3.8164 Angstrom in Bohr
    R0 = 7.21195
    # 1.0061 kJ/mol to Hartree
    epsilon = 0.000383203368

    # Avoid division by zero for gamma = 0.
    if gamma == 0.:
        alpha = 0.
    else:
        alpha = gamma / ((2**(-1/6) - (1 + (1 + gamma/epsilon)**0.5)**(-1/6)) * R0)
    
    def afir_func(coords3d):
        diffs = anp.diff(coords3d[inds], axis=1).reshape(-1, 3)
        rs = anp.linalg.norm(diffs, axis=1)

        omegas = (cov_rad_sums / rs)**p

        f = alpha * rho * (omegas*rs).sum() / omegas.sum()
        return f
    return afir_func


class AFIR(Calculator):

    def __init__(self, calculator, atoms, fragment_indices, gamma):
        super().__init__()

        self.calculator = calculator
        self.fragment_indices = fragment_indices
        self.gamma = gamma

        self.cov_radii = np.array([COVALENT_RADII[atom.lower()] for atom in atoms]) 
        self.afir_func = afir_closure(self.fragment_indices,
                                      self.cov_radii,
                                      self.gamma)
        self.afir_grad_func = autograd.grad(self.afir_func)

    def get_energy(self, atoms, coords):
        energy = self.calculator.get_energy(atoms, coords)["energy"]
        afir_energy = self.afir_func(coords.reshape(-1, 3))
        return {"energy": energy+afir_energy}

    def get_forces(self, atoms, coords):
        coords3d = coords.reshape(-1, 3)
        results = self.calculator.get_forces(atoms, coords)
        energy = results["energy"]
        forces = results["forces"]

        afir_energy = self.afir_func(coords3d)
        afir_forces = -self.afir_grad_func(coords3d).flatten()
        return {"energy": energy+afir_energy,
                "forces": forces+afir_forces,
        }
