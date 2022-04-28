import sys
import time

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.calculators import ORCA, XTB, ORCA5

CALC_CLASSES = {
    "orca": ORCA.ORCA,
    "orca5": ORCA5.ORCA5,
    "xtb": XTB.XTB,
}

try:
    from pysisyphus.calculators import PySCF

    CALC_CLASSES["pyscf"] = PySCF.PySCF
except ModuleNotFoundError:
    pass


def calcs_from_dict(calc_dict, base_name, calc_number, charge, mult, pal, mem):
    keys_calcs = {}
    calc_kwargs_ = {
        "calc_number": calc_number,
        "charge": charge,
        "mult": mult,
        "pal": pal,
        "mem": mem,
    }
    # Try to distinguish the calculators according to their base_name. If no
    # base_name is given we rely on the calc_number.
    base_name = base_name if base_name != "" else f"{calc_number:03d}"
    for i, (key, kwargs) in enumerate(calc_dict.items()):
        calc_kwargs = calc_kwargs_.copy()
        calc_kwargs["base_name"] = f"{base_name}_{key}"
        # Don't modify original dict
        kwargs = kwargs.copy()
        type_ = kwargs.pop("type")
        calc_kwargs.update(**kwargs)
        calc_cls = CALC_CLASSES[type_]
        keys_calcs[key] = calc_cls(**calc_kwargs)
    return keys_calcs


class MultiCalc(Calculator):
    def __init__(self, calcs, **kwargs):
        super().__init__(**kwargs)

        self.do_hess = {
            key: calc_dict.pop("do_hess", False) for key, calc_dict in calcs.items()
        }
        self.keys_calcs = calcs_from_dict(
            calcs,
            self.base_name,
            self.calc_number,
            self.charge,
            self.mult,
            self.pal,
            self.mem,
        )
        max_len = max([len(key) for key in self.keys_calcs.keys()])
        self.max_fmt = f" >{max_len+2}"

    def run_calculation(self, atoms, coords, **prepare_kwargs):
        all_results = {}
        for name, calc in self.keys_calcs.items():
            if self.do_hess[name]:
                run_func = calc.get_hessian
            else:
                run_func = calc.run_calculation
            key = f"{self.base_name}_{name}"
            print(f"Running {name:{self.max_fmt}} ... ", end="")
            try:
                start = time.time()
                results = run_func(atoms, coords, **prepare_kwargs)
                end = time.time()
                duration = end - start
                print(f"took {duration:.0f} s.")
            except Exception as err:
                print("\nCalculation failed!")
                print(err)
            sys.stdout.flush()
            all_results[key] = results
        return all_results
