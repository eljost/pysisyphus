from math import pi
import numpy as np

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import KB, AU2J


class LogFermi:
    def __init__(self, beta, radius, T=300, origin=(0.0, 0.0, 0.0)):
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
        t2 = np.exp(self.beta * (t1 - self.radius))

        energy = (self.kT * np.log(1 + t2)).sum()
        if not gradient:
            return energy

        gradient = self.kT * ((self.beta * t2) / ((1 + t2) * t1))[:, None] * t0
        return energy, gradient.flatten()

    def __repr__(self):
        return (
            f"LogFermi(beta={self.beta:.6f}, radius={self.radius:.6f}, "
            f"T={self.T:.6f}, origin={self.origin})"
        )


class HarmonicSphere:
    def __init__(self, k, radius, origin=(0.0, 0.0, 0.0)):
        self.k = k
        self.radius = radius
        self.origin = np.array(origin)

    def calc(self, coords3d, gradient=False):
        c3d_wrt_origin = coords3d - self.origin
        distances = np.linalg.norm(c3d_wrt_origin, axis=1)
        energies = np.where(distances > self.radius, self.k * distances ** 2, 0.0)
        energy = energies.sum()

        if not gradient:
            return energy

        """
        E(r(x)) = k*r**2
        dE(r(x))/dx = dE/dr * dr/dx
        dE/dr = 2*k*r
        dr/dx = x/r
        dE/dr * dr/dx = 2*k*x
        """
        gradient = np.where(distances > self.radius, 2 * self.k * c3d_wrt_origin, 0.0)

        return energy, gradient.flatten()

    @property
    def surface_area(self):
        """In Bohr**2"""
        return 4 * pi * self.radius**2

    def instant_pressure(self, coords3d):
        _, gradient = self.calc(coords3d, gradient=True)
        norm = np.linalg.norm(gradient)
        p = norm / self.surface_area
        return p


class ExternalPotential(Calculator):

    available_potentials = {
        "logfermi": LogFermi,
        "harmonic_sphere": HarmonicSphere,
    }

    def __init__(self, calculator=None, potentials=None, **kwargs):
        super().__init__(**kwargs)

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
        if self.calculator is not None:
            results = self.calculator.get_energy(atoms, coords)
        else:
            results = {"energy": 0.0}
        results["energy"] += potential_energy
        return results

    def get_potential_forces(self, coords):
        coords3d = coords.reshape(-1, 3)
        energies_gradients = [
            pot.calc(coords3d, gradient=True) for pot in self.potentials
        ]
        energies, gradients = zip(*energies_gradients)
        self.log(f"Energies from external potential: {energies}")
        energy = sum(energies)
        forces = -np.sum(gradients, axis=0)
        self.log(f"Forces from external potential: {forces}")
        return energy, forces

    def get_forces(self, atoms, coords):
        potential_energy, potential_forces = self.get_potential_forces(coords)
        if self.calculator is not None:
            results = self.calculator.get_forces(atoms, coords)
        else:
            results = {"energy": 0.0, "forces": np.zeros_like(coords)}
        results["energy"] += potential_energy
        results["forces"] += potential_forces
        return results

    def get_hessian(self, atoms, coords):
        raise Exception("Hessian is not implemented for ExternalPotential!")
