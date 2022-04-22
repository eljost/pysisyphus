from pysisyphus.constants import AU2KJPERMOL, AU2KCALPERMOL, AU2EV, ANG2BOHR, BOHR2M

import numpy as np
import yaml
from yaml.constructor import ConstructorError

_UNIT_MAP = {
    # Energies are converted to Hartree (E_h)
    "Eh": 1,
    "kJpermol": 1 / AU2KJPERMOL,  # kJ mol⁻¹
    "kcalpermol": 1 / AU2KCALPERMOL,  # kcal mol⁻¹
    "eV": 1 / AU2EV,  # eV
    # Times are converted to fs (1e-15 s)
    "fs": 1,
    "ps": 1e-3,
    "ns": 1e-6,
    # Lengths are converted to Bohr (a_0)
    "a0": 1,
    "angstrom": 1 / ANG2BOHR,
    "nm": 1 / BOHR2M * 1e9,
}
_UNITS = list(_UNIT_MAP.keys())
UNITS = _UNITS


def get_constructor(unit):
    conv = _UNIT_MAP[unit]

    def constructor(loader, node):
        try:
            data = float(loader.construct_scalar(node)) * conv
        except ConstructorError:
            data = np.array(loader.construct_sequence(node), dtype=float) * conv
        return data

    return constructor


def get_loader(units=_UNITS):
    loader = yaml.SafeLoader
    for unit in units:
        tag = f"!{unit}"
        loader.add_constructor(tag, get_constructor(unit))
    return loader
