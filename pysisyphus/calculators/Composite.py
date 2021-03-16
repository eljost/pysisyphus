from sympy import sympify

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.calculators import ORCA


CALC_CLASSES = {
    "orca": ORCA.ORCA,
}


class Composite(Calculator):
    def __init__(self, final, keys_calcs=None, from_dict=None, **kwargs):
        assert keys_calcs or from_dict
        super().__init__(**kwargs)

        if from_dict:
            keys_calcs = {}
            calc_kwargs_ = {
                "charge": self.charge,
                "mult": self.mult,
                "pal": self.pal,
                "mem": self.mem,
            }
            for i, (key, kwargs) in enumerate(from_dict.items()):
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

    def get_energy(self, atoms, coords, **prepare_kwargs):
        subst = self.final
        energies = {}
        for key, calc in self.keys_calcs.items():
            energy = calc.get_energy(atoms, coords, **prepare_kwargs)["energy"]
            energies[key] = energy

        final_energy = sympify(subst).subs(energies)
        results = {
            "energy": final_energy,
        }
        return results
