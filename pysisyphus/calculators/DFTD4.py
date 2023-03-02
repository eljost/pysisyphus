import dataclasses
from typing import Optional

try:
    from dftd4.interface import DampingParam, DispersionModel

    HAS_DFTD4 = True
except ModuleNotFoundError:
    HAS_DFTD4 = False
    # Dummy definition, so pysisyphus does not crash if another module tries to
    # import DFT4.
    DFTD4 = None

from pysisyphus.calculators.Calculator import Calculator
from pysisyphus.elem_data import nuc_charges_for_atoms


@dataclasses.dataclass
class D4Params:
    s6: float
    s8: float
    s9: float
    a1: float
    a2: float


OWN_METHODS = {
    # ORCA 5.0.3 ... your r²scan-3c bug ... why are you doing this to me? Why?
    "r2scan-3c": D4Params(s6=1.0, s8=0.0, s9=2.0, a1=0.42, a2=5.65),
}


class DFTD4(Calculator):
    def __init__(
        self,
        method: Optional[str] = None,
        damp_params: Optional[D4Params] = None,
        model_params: Optional[dict] = None,
        **kwargs
    ):
        assert HAS_DFTD4, "DFTD4 python package is not installed!"
        super().__init__(**kwargs)

        if model_params is None:
            model_params = dict()
        self.model_params = model_params

        if method in OWN_METHODS:
            damp_params = OWN_METHODS[method]

        if damp_params is not None:
            self.damp_param = DampingParam(**dataclasses.asdict(damp_params))
        # if method is not None:
        elif method is not None:
            self.method = method
            self.damp_param = DampingParam(method=method)
        else:
            raise Exception("Please provide either 'method' or 'params'!")

    def get_model(self, atoms, coords):
        atomic_numbers = nuc_charges_for_atoms(atoms)
        model = DispersionModel(
            atomic_numbers, coords.reshape(-1, 3), **self.model_params
        )
        return model

    def get_dispersion(self, atoms, coords, grad):
        model = self.get_model(atoms, coords)
        res = model.get_dispersion(self.damp_param, grad=grad)
        results = {
            "energy": float(res.get("energy")),
        }
        if grad:
            results["forces"] = -res.get("gradient").flatten()
        return results

    def get_energy(self, atoms, coords, **prepare_kwargs):
        return self.get_dispersion(atoms, coords, grad=False)

    def get_forces(self, atoms, coords, **prepare_kwargs):
        return self.get_dispersion(atoms, coords, grad=True)

    # This calculator relies on a finite differences numerical Hessian, as provided
    # by the Calculator parent class.
