import numpy as np
try:
    from xtb.interface import Environment, Param, Calculator as XTBCalculator
    from xtb.libxtb import VERBOSITY_MINIMAL, VERBOSITY_FULL, VERBOSITY_MUTED
except ModuleNotFoundError:
    pass

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.elem_data import ATOMIC_NUMBERS


class PyXTB(Calculator):
    avail_params = {
        "0": Param.GFN0xTB,
        "1": Param.GFN1xTB,
        "2": Param.GFN2xTB,
        "ff": Param.GFNFF,
    }

    avail_verbosities = {
        0: VERBOSITY_MUTED,
        1: VERBOSITY_MINIMAL,
        2: VERBOSITY_FULL,
    }

    def __init__(self, *args, gfn=2, verbosity=0, **kwargs):
        super().__init__(*args, **kwargs)

        self.env = Environment()
        self.gfn = str(gfn)
        self.verbosity = self.avail_verbosities[verbosity]
        self.param = self.avail_params[self.gfn]
        self.uhf = self.mult - 1

    def get_calculator(self, atoms, coords):
        numbers = np.array([ATOMIC_NUMBERS[atom.lower()] for atom in atoms])
        calc = XTBCalculator(
            self.param, numbers, coords.reshape(-1, 3), charge=self.charge, uhf=self.uhf
        )
        calc.set_verbosity(self.verbosity)
        return calc

    def get_energy(self, atoms, coords, **prepare_kwargs):
        calc = self.get_calculator(atoms, coords)
        xtb_result = calc.singlepoint(copy=True)
        results = {
            "energy": xtb_result.get_energy(),
        }
        return results

    def get_forces(self, atoms, coords):
        calc = self.get_calculator(atoms, coords)
        xtb_result = calc.singlepoint(copy=True)
        results = {
            "energy": xtb_result.get_energy(),
            "forces": -xtb_result.get_gradient().flatten(),
        }
        return results
