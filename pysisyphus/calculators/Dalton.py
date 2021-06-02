try:
    import daltonproject as dp
except ModuleNotFoundError:
    print("daltonproject is not installed!")

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.constants import BOHR2ANG


class Dalton(Calculator):

    conf_key = "dalton"

    def __init__(self, basis, method="hf", **kwargs):
        super().__init__(**kwargs)

        self._basis = basis
        self._method = method

        self.basis = dp.Basis(basis=self._basis)
        self.method = dp.QCMethod(self._method)
        self.compute_settings = dp.ComputeSettings(mpi_num_procs=1)

    def prepare_input(self, atoms, coords):
        coords3d = coords.reshape(-1, 3) * BOHR2ANG
        dp_atoms = "; ".join(
            [f"{atom} {x} {y} {z}" for atom, (x, y, z) in zip(atoms, coords3d)]
        )
        molecule = dp.Molecule(atoms=dp_atoms, charge=self.charge)

        return molecule

    def compute(self, mol, prop):
        return dp.dalton.compute(
            mol, self.basis, self.method, prop, compute_settings=self.compute_settings
        )

    def get_energy(self, atoms, coords, **prepare_kwargs):
        mol = self.prepare_input(atoms, coords)
        prop = dp.Property(energy=True)
        res = self.compute(mol, prop)

        results = {
            "energy": res.energy,
        }
        return results


    def get_forces(self, atoms, coords, **prepare_kwargs):
        mol = self.prepare_input(atoms, coords)
        prop = dp.Property(energy=True, gradients=True)
        res = self.compute(mol, prop)

        results = {
            "energy": res.energy,
            "forces": -res.gradients.flatten(),
        }
        return results

    def __str__(self):
        return f"Dalton({self.name})"
