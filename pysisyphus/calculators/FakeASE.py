from pysisyphus.constants import BOHR2ANG

class FakeASE:

    def __init__(self, calc):
        self.calc = calc

        self.results = dict()

    def get_atoms_coords(self, atoms):
        return (atoms.get_chemical_symbols(),
                # Convert ASE Angstrom to Bohr for pysisyphus
                atoms.get_positions().flatten() / BOHR2ANG
        )

    def get_potential_energy(self, atoms=None):
        atoms, coords = self.get_atoms_coords(atoms)
        results = self.calc.get_energy(atoms, coords)

        energy = results["energy"]
        # self.results["energy"] = results["energy"]

        return energy

    def get_forces(self, atoms=None):
        atoms, coords = self.get_atoms_coords(atoms)
        results = self.calc.get_forces(atoms, coords)

        # Convert back to Angstrom for ASE
        forces = results["forces"].reshape(-1, 3) / BOHR2ANG

        # self.results["energy"] = results["energy"]
        # self.results["forces"] = forces

        return forces
