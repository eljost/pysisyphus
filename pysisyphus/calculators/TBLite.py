import numpy as np

try:
    from tblite.interface import Calculator as tbCalculator

    can_pyxtb = True
except ModuleNotFoundError:
    can_pyxtb = False

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.elem_data import ATOMIC_NUMBERS


class TBLite(Calculator):

    def __init__(
        self,
        *args,
        gfn: int = 2,
        acc: float = 0.1,
        verbosity: int = 0,
        keep_calculator: bool = False,
        **kwargs
    ):
        super().__init__(*args, **kwargs)

        self.gfn = str(gfn)
        self.acc = float(acc)
        self.verbosity = int(verbosity)
        self.keep_calculator = bool(keep_calculator)
        self._calculator = None

        self.param = {
            "1": "GFN1-xTB",
            "2": "GFN2-xTB",
        }
        self.method = self.param[self.gfn]
        self.uhf = self.mult - 1

    def get_calculator(self, atoms, coords):
        # Reuse old calculator and updatec coordinates
        if self._calculator:
            self._calculator.update(coords.reshape(-1, 3))
            return self._calculator

        # Create new calculator
        numbers = np.array([ATOMIC_NUMBERS[atom.lower()] for atom in atoms])
        calc = tbCalculator(
            self.method,
            numbers,
            coords.reshape(-1, 3),
            charge=self.charge,
            uhf=self.uhf,
        )
        calc.set("verbosity", self.verbosity)
        calc.set("accuracy", self.acc)

        # Keep calculator, if requested
        if self.keep_calculator and (self._calculator is None):
            self._calculator = calc
        return calc

    def get_energy(self, atoms, coords, **prepare_kwargs):
        calc = self.get_calculator(atoms, coords)
        xtb_result = calc.singlepoint(copy=True)
        results = {
            "energy": xtb_result.get("energy"),
        }
        return results

    def get_forces(self, atoms, coords, **prepare_kwargs):
        calc = self.get_calculator(atoms, coords)
        xtb_result = calc.singlepoint(copy=True)
        results = {
            "energy": xtb_result.get("energy"),
            "forces": -xtb_result.get("gradient").flatten(),
        }
        return results
