import numpy as np

try:
    from simtk.openmm import app
    from simtk import openmm as mm
    from simtk import unit
except ModuleNotFoundError:
    print("OpenMM is not installed!")

from pysisyphus.constants import BOHR2ANG


class OpenMM:
    def __init__(self, topology, params):
        self.topology = topology
        self.params = params
        self.system = self.topology.createSystem(self.params, nonbondedMethod=app.NoCutoff)

        dummy_int = mm.VerletIntegrator(0.001)
        self.simulation = app.Simulation(self.topology, self.system, dummy_int)
        self.context = self.simulation.context

        self.en2au = 1 / unit.hartree / unit.AVOGADRO_CONSTANT_NA
        self.frc2au = 1 / (unit.hartree / unit.bohr) / unit.AVOGADRO_CONSTANT_NA

    def get_energy(self, atoms, coords, **prepare_kwargs):
        results = self.get_forces(atoms, coords, **prepare_kwargs)
        return {"energy": results["energy"]}

    def get_forces(self, atoms, coords, **prepare_kwargs):
        positions = coords.reshape(-1, 3) * BOHR2ANG / 10  # 2d array in nm
        self.context.setPositions(positions)
        state = self.context.getState(getEnergy=True, getForces=True)
        energy = state.getPotentialEnergy() * self.en2au
        forces = state.getForces(asNumpy=True)
        forces = (forces * self.frc2au).reshape(-1)

        energy = float(energy)
        forces = np.array(forces)
        return {"energy": energy, "forces": forces}
