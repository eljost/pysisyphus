import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import KB, AU2J


class LogFermi:

    def __init__(self, beta, radius, T=300, origin=(0., 0., 0.)):
        """As described in the XTB docs.

        https://xtb-docs.readthedocs.io/en/latest/xcontrol.html#confining-in-a-cavity
        """
        self.beta = beta
        self.radius = radius
        self.T = T
        self.origin = np.array(origin)

        # In Hartree
        self.kT = KB * self.T / AU2J

    def calc(self, coords3d, gradient=False):
        t0 = coords3d - self.origin[None, :]
        t1 = np.linalg.norm(t0, axis=1)
        t2 = np.exp(self.beta*(t1 - self.radius))

        energy = (self.kT * np.log(1 + t2)).sum()
        if not gradient:
            return energy

        gradient = self.kT * ((self.beta*t2) / ((1+t2) * t1))[:, None] * t0
        return energy, gradient.flatten()

    def __repr__(self):
        return f"LogFermi(beta={self.beta:.6f}, radius={self.radius:.6f}, " \
               f"T={self.T:.6f}, origin={self.origin})"


class ExternalPotential(Calculator):

    available_potentials = {
        "logfermi": LogFermi,
    }

    def __init__(self, calculator, *args, potentials=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.calculator = calculator

        self.potentials = list()
        self.log("Creating external potentials")
        for i, pot_kwargs in enumerate(potentials):
            pot_key = pot_kwargs.pop("type")
            pot_cls = self.available_potentials[pot_key]
            pot = pot_cls(**pot_kwargs)
            self.potentials.append(pot)
            self.log(f"\t{i:02d}: {pot}")

    def get_potential_energy(self, coords):
        coords3d = coords.reshape(-1, 3)
        potential_energies = [pot.calc(coords3d) for pot in self.potentials]
        potential_energy = sum(potential_energies)
        self.log(f"Energies from external potential: {potential_energies}")
        return potential_energy

    def get_energy(self, atoms, coords):
        potential_energy = self.get_potential_energy(coords)
        results = self.calculator.get_energy(atoms, coords)
        results["energy"] += potential_energy
        return results

    def get_potential_forces(self, coords):
        coords3d = coords.reshape(-1, 3)
        energies_gradients = [pot.calc(coords3d, gradient=True)
                              for pot in self.potentials]
        energies, gradients = zip(*energies_gradients)
        self.log(f"Energies from external potential: {energies}")
        energy = sum(energies)
        forces = -np.sum(gradients, axis=0)
        self.log(f"Forces from external potential: {forces}")
        return energy, forces

    def get_forces(self, atoms, coords):
        potential_energy, potential_forces = self.get_potential_forces(coords)
        results = self.calculator.get_forces(atoms, coords)
        results["energy"] += potential_energy
        results["forces"] += potential_forces
        return results

    def get_hessian(self, atoms, coords):
        raise Exception("Hessian is not implemented for ExternalPotential!")
