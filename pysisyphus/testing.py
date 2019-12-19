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
            cmd = Config[calculator]["cmd"]
            # _ = shutil.which(cmd)
            # print("\t", calculator, _)
            available = bool(shutil.which(cmd))
        except KeyError:
            pass

        if calculator == "pyscf":
            try:
                import pyscf
                available = True
            except ImportError:
                pass

        reason = _reason.format(calculator)
        _using_cache[calculator] = pytest.mark.skipif(not available, reason=reason)
    return _using_cache[calculator]



using_pytest = pytest.mark.skipif(not using("pyscf"), reason=_reason.format("pyscf"))
