import os
import time
import timeit

import numpy as np
import pytest

from pysisyphus.wavefunction import Wavefunction
from pysisyphus.testing import using


@pytest.mark.skip
@using("numba")
def test_integral_times():
    from pysisyphus.wavefunction.backend_numba import (
        _NUMBA_MODULES,
        get_func_data,
    )

    os.environ["NUMBA_DEBUG_CACHE"] = "1"
    print("Timing numba integral imports/jitting")
    for key in _NUMBA_MODULES.keys():
        print(f"Importing {key} ... ", end="")
        dur = timeit.timeit(lambda: get_func_data(key))
        print(f"this took {dur: >12.4f} s")
    try:
        del os.environ["NUMBA_DEBUG_CACHE"]
    except KeyError:
        pass


@pytest.fixture(scope="module")
def numba_wf():
    return Wavefunction.from_file(
        "lib:01_ch4_tzvp.json",
        shell_kwargs={
            "backend": "numba",
        },
    )


@pytest.fixture
def py_wf():
    return Wavefunction.from_file("lib:01_ch4_tzvp.json")


@using("numba")
@pytest.mark.parametrize(
    "func, kwargs",
    (
        ("get_T_cart", None),
        ("get_dipole_ints_cart", {"origin": np.zeros(3)}),
        ("get_diag_quadrupole_ints_cart", {"origin": np.zeros(3)}),
        ("get_quadrupole_ints_cart", {"origin": np.zeros(3)}),
    ),
)
def test_numba_integral_matrices(func, kwargs, py_wf, numba_wf):
    if kwargs is None:
        kwargs = {}
    py_shells = py_wf.shells

    py_dur = time.time()
    py_ints = getattr(py_shells, func)(**kwargs)
    py_dur = time.time() - py_dur

    numba_shells = numba_wf.shells

    # Numba warump
    getattr(numba_shells, func)(**kwargs)

    numba_dur = time.time()
    numba_ints = getattr(numba_shells, func)(**kwargs)
    numba_dur = time.time() - numba_dur

    print(f"python: {py_dur:.4f} s, numba: {numba_dur:.4f} s")
    np.testing.assert_allclose(numba_ints, py_ints, atol=1e-14)
