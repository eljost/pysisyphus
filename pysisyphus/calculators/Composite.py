import numpy as np
import sympy as sym

# from sympy import sympify, lambdify,

from pysisyphus.calculators.Calculator import Calculator

from pysisyphus.calculators import ORCA


CALC_CLASSES = {
    "orca": ORCA.ORCA,
}
try:
    from pysisyphus.calculators import PySCF

    CALC_CLASSES["pyscf"] = PySCF.PySCF
except ModuleNotFoundError:
    pass


class Composite(Calculator):
    def __init__(self, final, keys_calcs=None, calcs=None, remove_translation=False,
            **kwargs):
        # Either directly supply a dictionary with Calculator objects (key_calcs)
        # or a dictionary containing information to set up calculators.
        assert keys_calcs or calcs

        super().__init__(**kwargs)

        if calcs:
            keys_calcs = {}
            calc_kwargs_ = {
                "charge": self.charge,
                "mult": self.mult,
                "pal": self.pal,
                "mem": self.mem,
            }
            for i, (key, kwargs) in enumerate(calcs.items()):
                calc_kwargs = calc_kwargs_.copy()
                calc_kwargs["calc_number"] = i
                calc_kwargs["base_name"] = f"composite_{key}"
                # Don't modify original dict
                kwargs = kwargs.copy()
                type_ = kwargs.pop("type")
                calc_kwargs.update(**kwargs)
                calc_cls = CALC_CLASSES[type_]
                keys_calcs[key] = calc_cls(**calc_kwargs)

        self.keys_calcs = keys_calcs
        assert all([key in final for key in self.keys_calcs.keys()])
        self.final = final
        self.remove_translation = remove_translation

        self.energy_expr = sym.sympify(self.final)
        self.forces_args = sym.symbols(" ".join(self.keys_calcs.keys()))
        self.forces_expr = sym.lambdify(self.forces_args, self.energy_expr)

    def get_energy(self, atoms, coords, **prepare_kwargs):
        energies = {}
        for key, calc in self.keys_calcs.items():
            energy = calc.get_energy(atoms, coords, **prepare_kwargs)["energy"]
            energies[key] = energy

        final_energy = self.energy_expr.subs(energies).evalf()
        results = {
            "energy": final_energy,
        }
        return results

    def get_forces(self, atoms, coords, **prepare_kwargs):
        energies = {}
        forces = {}
        for key, calc in self.keys_calcs.items():
            results = calc.get_forces(atoms, coords, **prepare_kwargs)
            energies[key] = results["energy"]
            forces[key] = results["forces"]
        keys = self.keys_calcs.keys()
        for key in keys:
            self.log(f"|forces_{key}|={np.linalg.norm(forces[key]):.6f}")
        self.log("")

        final_energy = self.energy_expr.subs(energies).evalf()
        final_forces = self.forces_expr(**forces)

        # Remove overall translation
        if self.remove_translation:
            f3d = final_forces.reshape(-1, 3)
            f3d -= f3d.mean(axis=0)[None, :]

        results = {
            "energy": final_energy,
            "forces": final_forces,
        }
        return results

    def run_calculation(self, atoms, coords, **prepare_kwargs):
        return self.get_energy(atoms, coords, **prepare_kwargs)
