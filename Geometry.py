class Geometry:

    def __init__(self, atoms, coords):
        self.atoms = atoms
        self.coords = coords

        self._energy = None
        self._forces = None
        self._hessian = None

    def set_calculator(self, calculator):
        self.calculator = calculator

    @property
    def energy(self):
        if not self._energy:
            results = self.calculator.get_energy(self.coords)
            self.set_results(results)
        return self._energy

    @energy.setter
    def energy(self, energy):
        self._energy = energy

    @property
    def forces(self):
        if not self._forces:
            results = self.calculator.get_forces(self.coords)
            self.set_results(forces)
        return self._forces

    @forces.setter
    def forces(self, forces):
        self._forces = forces

    @property
    def hessian(self):
        if not self._hessian:
            results = self.calculator.get_hessian(self.coords)
            self.set_results(results)
        return self._hessian

    @hessian.setter
    def hessian(self, hessian):
        self._hessian = hessian

    def set_results(self, results):
        for key in results:
            setattr(self, key, results[key])

    def calc_energy_and_forces(self):
        results = self.calculator.get_forces(self.coords)
        self.set_results(results)

    def __str__(self):
        return "{} atoms".format(len(self.atoms))
