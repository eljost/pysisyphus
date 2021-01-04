import importlib

import pytest
import shutil

from pysisyphus.config import Config, DEFAULTS


"""Inspired by
https://github.com/MolSSI/QCEngine/blob/master/qcengine/testing.py

Maybe implement something like this to make test access easier
    https://stackoverflow.com/a/36601576/12486216
"""


_reason = "Calculator {} is not available."
_using_cache = dict()
# Python modules to be imported using importlib.import_module
IMPORT_DICT = {
    "pyscf": "pyscf",
    "dalton": "daltonproject",
    "qcengine": "qcengine",
    "openmm": "simtk.openmm",
    "thermoanalysis": "thermoanalysis",
    "pyxtb": "xtb.interface",
    "obabel": "openbabel.openbabel",
}


def module_available(calculator):
    try:
        import_name = IMPORT_DICT[calculator]
    except KeyError:
        return False

    try:
        _ = importlib.import_module(import_name)
        available = True
    except (ModuleNotFoundError, ImportError):
        available = False
    return available


def using(calculator):
    calculator = calculator.lower()

    if calculator not in _using_cache:
        available = False

        # Look into .pysisyphusrc first
        try:
            cmd = Config[calculator]["cmd"]
        except KeyError:
            # Try defaults last
            try:
                cmd = DEFAULTS[calculator]
            except KeyError:
                cmd = None

            # Look for dscf availability
            if calculator == "turbomole":
                cmd = "dscf"

        # Check if cmd is available on $PATH
        if cmd is not None:
            available = bool(shutil.which(cmd))
        # Handling native python packages from here
        else:
            available = module_available(calculator)

        reason = _reason.format(calculator)
        _using_cache[calculator] = pytest.mark.skipif(not available, reason=reason)
    return _using_cache[calculator]


def available(calculator):
    # True when skipif is False
    return not using(calculator).args[0]


using_pyscf = using("pyscf")
using_gaussian16 = using("gaussian16")
using_turbomole = using("turbomole")
using_qcengine = using("qcengine")
