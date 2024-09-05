import importlib

import numba
from numba import i8, f8
from numba.experimental import jitclass
from numba.core import types
from numba.core.extending import overload_method
from numba.experimental import structref
import numpy as np


shell_spec = [
    ("L", i8),
    ("center", f8[:]),
    ("coeffs", f8[:]),
    ("exps", f8[:]),
    ("index", i8),
    ("size", i8),
]


@jitclass(shell_spec)
class NumbaShell(object):
    def __init__(self, L, center, coeffs, exps, index, size):
        self.L = L
        self.center = center
        self.coeffs = coeffs
        self.exps = exps
        self.index = index
        self.size = size

    def as_tuple(self):
        return self.L, self.center, self.coeffs, self.exps, self.index, self.size

    @property
    def cart_size(self):
        return (self.L + 2) * (self.L + 1) // 2

    @property
    def sph_size(self):
        return 2 * self.L + 1


_NUMBA_MODULES = {
    "int1e_ovlp": ("ovlp3d", 0),
    "int1e_r": ("dipole3d", 3),
    "int1e_rr": ("quadrupole3d", 6),
    "int1e_drr": ("diag_quadrupole3d", 3),
    "int1e_kin": ("kinetic3d", 0),
}

# This dict will be populated as needed, as one import costs about ~2s
_func_data = {}


def get_func_data(key):
    """Get func_dict and component.

    _func_data dict is laziliy populated, as there is still a ~ 2 s compile time
    per import :/
    """
    if key not in _func_data:
        module_name, components = _NUMBA_MODULES[key]
        full_module_name = f"pysisyphus.wavefunction.ints_numba.{module_name}"
        module = importlib.import_module(full_module_name)
        func_dict = getattr(module, "get_func_dict")()
        _func_data[key] = (func_dict, components)
    return _func_data[key]


_R0 = np.zeros(3)


@numba.jit(parallel=True, nopython=True, cache=True)
def get_2c_ints_cart(
    shells_a,
    shells_b,
    func_dict,
    components,
    symmetric=True,
    R=_R0,
):
    tot_size_a = 0
    for shell in shells_a:
        tot_size_a += shell.size

    tot_size_b = 0
    for shell in shells_b:
        tot_size_b += shell.size

    # Allocate final integral array
    integrals = np.zeros((components, tot_size_a, tot_size_b))
    shells_b = shells_a
    nshells_a = len(shells_a)
    nshells_b = len(shells_b)

    # Start loop over contracted gaussians in shells_a
    for i in numba.prange(nshells_a):
        shell_a = shells_a[i]
        La, A, das, axs, indexa, sizea = shell_a.as_tuple()
        na = len(das)
        slicea = slice(indexa, indexa + sizea)
        # Start loop over contracted gaussians in shells_b
        for j in range(i, nshells_b):
            shell_b = shells_b[j]
            Lb, B, dbs, bxs, indexb, sizeb = shell_b.as_tuple()
            nb = len(dbs)
            sliceb = slice(indexb, indexb + sizeb)
            result = np.zeros(max(components, 1) * sizea * sizeb)
            # Pick correct function depending on La and Lb
            func = func_dict[(La, Lb)]

            # Start loop over primitives
            for k in range(na):
                ax = axs[k]
                da = das[k]
                for l in range(nb):
                    bx = bxs[l]
                    db = dbs[l]
                    func(ax, da, A, bx, db, B, R, result)
            integrals[:, slicea, sliceb] = result.reshape(components, sizea, sizeb)
            # End loop over primitives gaussians

            if symmetric and (i != j):
                for k in range(indexa, indexa + sizea):
                    for l in range(indexb, indexb + sizeb):
                        integrals[:, l, k] = integrals[:, k, l]
        # End loop over contracted gaussians in shells_b
    # End loop over contracted gaussians in shells_a
    return integrals


def to_numba_shells(shells):
    numba_shells = numba.typed.List()
    for shell in shells:
        numba_shell = NumbaShell(*shell.as_tuple())
        numba_shells.append(numba_shell)
    return numba_shells


"""
@numba.jit(nopython=True, cache=True)
def to_numba_shells_from_tuples(shells):
    numba_shells = numba.typed.List()
    for shell in shells:
        numba_shell = NumbaShell(*shell)
        numba_shells.append(numba_shell)
    return numba_shells
"""


@structref.register
class ShellStructType(types.StructRef):
    def preprocess_fields(self, fields):
        # This method is called by the type constructor for additional
        # preprocessing on the fields.
        # Here, we don't want the struct to take Literal types.
        return tuple((name, types.unliteral(typ)) for name, typ in fields)


SST = ShellStructType(
    [
        ("L", numba.i8),
        ("center", numba.f8[:]),
        ("center_ind", numba.i8),
        ("coeffs", numba.f8[:]),
        ("exps", numba.f8[:]),
        ("index", numba.i8),
        ("size", numba.i8),
    ]
)


class ShellStruct(structref.StructRefProxy):
    def __new__(cls, L, center, center_ind, coeffs, exps, index, size):
        return structref.StructRefProxy.__new__(
            cls, L, center, center_ind, coeffs, exps, index, size
        )


@overload_method(ShellStructType, "cart_size")
def ol_cart_size(self):
    def inner(self):
        return (self.L + 2) * (self.L + 1) // 2

    return inner


@overload_method(ShellStructType, "sph_size")
def ol_sph_size(self):
    def inner(self):
        return 2 * self.L + 1

    return inner


@overload_method(ShellStructType, "as_tuple")
def ol_as_tuple(self):
    def inner(self):
        return (
            self.L,
            self.center,
            self.center_ind,
            self.coeffs,
            self.exps,
            self.index,
            self.size,
        )

    return inner


structref.define_proxy(
    ShellStruct,
    ShellStructType,
    ["L", "center", "center_ind", "coeffs", "exps", "index", "size"],
)


def to_numba_shellstructs(shells):
    shellstructs = numba.typed.List()
    for shell in shells:
        L, center, coeffs, exps, index, size = shell.as_tuple()
        center_ind = shell.center_ind
        mst = ShellStruct(L, center, center_ind, coeffs, exps, index, size)
        shellstructs.append(mst)
    return shellstructs


def get_1el_ints_cart(
    shells_a,
    func_dict,
    shells_b=None,
    components=0,
    R=np.zeros(3),
    **kwargs,
):
    org_components = components
    components = max(1, components)
    symmetric = shells_b is None
    if symmetric:
        shells_b = shells_a

    integrals = get_2c_ints_cart(
        shells_a,
        shells_b,
        func_dict,
        components=components,
        symmetric=symmetric,
        R=R,
    )
    if org_components == 0:
        integrals = np.squeeze(integrals, axis=0)
    return integrals
