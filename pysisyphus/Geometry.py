import logging

import numpy as np

from pysisyphus.constants import BOHR2ANG
from pysisyphus.xyzloader import make_xyz_str
from pysisyphus.masses import MASS_DICT

class Geometry:

    def __init__(self, atoms, coords):
        self.atoms = atoms
        self._coords = coords

        self._energy = None
        self._forces = None
        self._hessian = None

        self.masses = [MASS_DICT[atom.lower()] for atom in self.atoms]
        # Some of the analytical potentials are only 2D
        repeat_masses = 2 if (self._coords.size == 2) else 3
        self.masses_rep = np.repeat(self.masses, repeat_masses)

    def clear(self):
        self.calculator = None
        self._energy = None
        self._forces = None
        self._hessian = None

    def set_calculator(self, calculator):
        self.clear()
        self.calculator = calculator

    @property
    def mm_inv(self):
        return np.linalg.inv(np.diag(self.masses_rep))

    @property
    def mm_sqrt_inv(self):
        return np.linalg.inv(np.sqrt(np.diag(self.masses_rep)))

    @property
    def coords(self):
        return self._coords

    @coords.setter
    def coords(self, coords):
        self._coords = coords
        # Reset all values because no calculations with the new coords
        # have been performed yet.
        self._energy = None
        self._forces = None
        self._hessian = None

    @property
    def mw_coords(self):
        return np.sqrt(self.masses_rep) * self._coords

    @mw_coords.setter
    def mw_coords(self, mw_coords):
        self.coords = mw_coords / np.sqrt(self.masses_rep)

    @property
    def energy(self):
        if self._energy is None:
            results = self.calculator.get_energy(self.atoms, self.coords)
            self.set_results(results)
        return self._energy

    @energy.setter
    def energy(self, energy):
        self._energy = energy

    @property
    def forces(self):
        if self._forces is None:
            results = self.calculator.get_forces(self.atoms, self.coords)
            self.set_results(results)
        return self._forces

    @forces.setter
    def forces(self, forces):
        self._forces = forces

    @property
    def gradient(self):
        return -self.forces

    @gradient.setter
    def gradient(self, gradient):
        self._forces = -gradient

    @property
    def mw_gradient(self):
        return -self.forces / np.sqrt(self.masses_rep)

    @property
    def hessian(self):
        if self._hessian is None:
            results = self.calculator.get_hessian(self.atoms, self.coords)
            self.set_results(results)
        return self._hessian

    @property
    def mw_hessian(self):
        # M^(-1/2) H M^(-1/2)
        return self.mm_sqrt_inv.dot(self.hessian).dot(self.mm_sqrt_inv)

    @hessian.setter
    def hessian(self, hessian):
        self._hessian = hessian

    def calc_energy_and_forces(self):
        results = self.calculator.get_forces(self.atoms, self.coords)
        self.set_results(results)

    def set_results(self, results):
        for key in results:
            setattr(self, key, results[key])
        self.results = results

    def as_xyz(self, comment=""):
        coords = self.coords * BOHR2ANG
        return make_xyz_str(self.atoms, coords.reshape((-1,3)), comment)

    def __str__(self):
        return "{} atoms".format(len(self.atoms))
