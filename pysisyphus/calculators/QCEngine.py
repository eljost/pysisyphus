import numpy as np
import qcengine as qcng
import qcelemental as qcel

from pysisyphus.calculators.Calculator import Calculator


class QCEngine(Calculator):

    def __init__(self, program, model, keywords=None, connectivity=None,
                 bond_order=1, **kwargs):
        super().__init__(**kwargs)

        self.program = program
        self.model = model
        if keywords is None:
            keywords = dict()
        self.keywords = dict(keywords)
        self.connectivity = connectivity
        self.bond_order = int(bond_order)

        # TODO: pal, memory

    def get_molecule(self, atoms, coords):
        c3d = coords.reshape(-1, 3)# * BOHR2ANG
        mol_inp = {
            "symbols": atoms,
            "geometry": c3d,
            "molecular_charge": self.charge,
            "molecular_multiplicity": self.mult,
        }

        if self.program == "openmm":
            if self.connectivity is None:
                self.log( "No connectivity specified! Using hacky connectivity "
                         f"guess with bond-order={self.bond_order}")
                connectivity = qcel.molutil.guess_connectivity(atoms, c3d,
                                                               threshold=1.3)
                self.connectivity = [(at1, at2, self.bond_order)
                                     for at1, at2 in connectivity]
            mol_inp["connectivity"] = self.connectivity

        molecule = qcel.models.Molecule.from_data(mol_inp)

        return molecule

    def compute(self, molecule, driver):
        inp = {
            "molecule": molecule,
            "driver": driver,
            "model": self.model,
            "keywords": self.keywords,
        }

        result = qcng.compute(inp, program=self.program, return_dict=True)
        self.calc_counter += 1
        return result

    def keep_stdout(self, res):
        if "stdout" in res:
            with open(self.make_fn(f"qce_{self.program}_stdout"), "w") as handle:
                    handle.write(res["stdout"])

    def get_energy(self, atoms, coords, **prepare_kwargs):
        mol = self.get_molecule(atoms, coords)

        res = self.compute(mol, driver="energy")
        self.keep_stdout(res)

        results = {
            "energy": res["return_result"],
        }

        return results

    def get_forces(self, atoms, coords, **prepare_kwargs):
        mol = self.get_molecule(atoms, coords)

        res = self.compute(mol, driver="gradient")
        self.keep_stdout(res)

        results = {
            "energy": res["properties"]["return_energy"],
            "forces": -np.array(res["return_result"]),
        }

        return results

    def get_hessian(self, atoms, coords, **prepare_kwargs):
        mol = self.get_molecule(atoms, coords)

        res = self.compute(mol, driver="hessian")
        self.keep_stdout(res)

        size = 3 * len(atoms)
        shape = (size, size)

        results = {
            "energy": res["properties"]["return_energy"],
            "hessian": np.array(res["return_result"]).reshape(shape),
        }

        return results

    def __str__(self):
        return f"QCECalculator({self.name})"
