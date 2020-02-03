import pytest
import shutil

from pysisyphus.config import Config


"""Inspired by
https://github.com/MolSSI/QCEngine/blob/master/qcengine/testing.py

Maybe implement something like this to make test access easier
    https://stackoverflow.com/a/36601576/12486216
"""


_reason = "Calculator {} is not available."
_using_cache = dict()


def using(calculator):
    calculator = calculator.lower()

    if calculator not in _using_cache:
        available = False

        try:
            # Look for dscf availability
            if calculator == "turbomole":
                cmd = "dscf"
            else:
                # Try to read from .pysisyphusrc
                cmd = Config[calculator]["cmd"]


            available = bool(shutil.which(cmd))
        except KeyError:
            pass

        if calculator == "pyscf":
            try:
                import pyscf
                available = True
            except ImportError:
                pass

        if calculator == "qcengine":
            try:
                import qcengine
                available = True
            except ImportError:
                pass

        if calculator == "openmm":
            try:
                import simtk.openmm
                available = True
            except ImportError:
                pass

        reason = _reason.format(calculator)
        _using_cache[calculator] = pytest.mark.skipif(not available, reason=reason)
    return _using_cache[calculator]



using_pyscf = using("pyscf")
using_gaussian16 = using("gaussian16")
using_turbomole = using("turbomole")
using_qcengine = using("qcengine")
