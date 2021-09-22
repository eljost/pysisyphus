import numpy as np

try:
    from simtk.openmm import app
    from simtk import openmm as mm
    from simtk import unit
except ModuleNotFoundError:
    print("OpenMM is not installed!")

from pysisyphus.constants import BOHR2ANG
from pysisyphus.calculators.Calculator import Calculator


class OpenMM(Calculator):
    def __init__(self, topology, params, **kwargs):
        super().__init__(**kwargs)

        self.topology = topology
        self.params = params
        self.system = self.topology.createSystem(
            self.params, nonbondedMethod=app.NoCutoff
        )

        dummy_int = mm.VerletIntegrator(0.001)
        self.simulation = app.Simulation(self.topology, self.system, dummy_int)
        self.context = self.simulation.context

        self.en2au = 1 / unit.hartree / unit.AVOGADRO_CONSTANT_NA
        self.frc2au = 1 / (unit.hartree / unit.bohr) / unit.AVOGADRO_CONSTANT_NA

    def get_charges(self):
        (nb_force,) = [
            force
            for force in self.system.getForces()
            if force.__class__.__name__ == "NonbondedForce"
        ]
        charges = list()
        for i in range(nb_force.getNumParticles()):
            charge, *_ = nb_force.getParticleParameters(i)
            charges.append(charge.value_in_unit(unit.elementary_charge))
        charges = np.array(charges)
        return charges

    def report_charges(self):
        """Total charge of the system.

        Adapted from modeller.py _addIons() in OpenMM.
        """

        charges = self.get_charges()
        total_charge = charges.sum()
        charges_str = np.array2string(charges, precision=4)
        self.log(f"OpenMM charges: {charges_str}")
        self.log(f"OpenMM system total charge: {total_charge:.6f} e")

    def get_dipole_moment(self, coords3d, reference=None, masses=None):
        def get_com(coords3d):
            total_mass = sum(masses)
            com = (1 / total_mass * coords3d * masses[:, None]).sum(axis=0)
            return com

        charges = self.get_charges()
        ref_dict = {
            None: lambda c3d: c3d,
            "centroid": lambda c3d: c3d.mean(axis=1)[None, :],
            "com": lambda c3d: c3d - get_com(c3d)[None, :],
        }
        coords3d = ref_dict[reference](coords3d)
        # Dipole moment vector in bohr * elementary_charge
        dip_moms = (coords3d * charges[:, None]).sum(axis=0)
        tot_dip_mom = np.linalg.norm(dip_moms)
        return dip_moms, tot_dip_mom

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
