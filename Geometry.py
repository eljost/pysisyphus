class Geometry:

    def __init__(self, coords):
        self.coords = coords

        self._energy = None
        self._forces = None
        self._hessian = None

    def set_calculator(self, calculator):
        self.calculator = calculator

    @property
    def energy(self):
        if not self._energy:
            self._energy = self.calculator.get_energy(self.coords)
        return self._energy

    @property
    def forces(self):
        if not self._forces:
            self._forces = self.calculator.get_forces(self.coords)
        return self._forces

    @property
    def hessian(self):
        if not self._hessian:
            self._hessian = self.calculator.get_hessian(self.coords)
        return self._hessian
