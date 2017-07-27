import numpy as np

class Geometry:

    def __init__(self, atoms, coords):
        self.atoms = atoms
        self._coords = coords

        self._energy = None
        self._forces = None
        self._hessian = None

    def set_calculator(self, calculator):
        self.calculator = calculator

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
    def hessian(self):
        if self._hessian is None:
            results = self.calculator.get_hessian(self.coords)
            self.set_results(results)
        return self._hessian

    @hessian.setter
    def hessian(self, hessian):
        self._hessian = hessian


    def calc_energy_and_forces(self):
        results = self.calculator.get_forces(self.atoms, self.coords)
        self.set_results(results)

    def set_results(self, results):
        for key in results:
            setattr(self, key, results[key])

    def __str__(self):
        return "{} atoms".format(len(self.atoms))
