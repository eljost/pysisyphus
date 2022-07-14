import json

import numpy as np

from pysisyphus.config import BASIS_LIB_DIR
from pysisyphus.elem_data import ATOMIC_NUMBERS
from pysisyphus.wavefunction import Shell, Shells


def basis_from_file(name):
    basis_path = BASIS_LIB_DIR / f"{name}.json"
    with open(basis_path) as handle:
        data = json.load(handle)
    elements = data["elements"]
    return elements


def shells_with_basis(atoms, coords, basis=None, name=None, **kwargs):
    assert (basis is not None) or (name is not None)
    if name is not None:
        basis = basis_from_file(name)

    coords3d = np.reshape(coords, (len(atoms), 3))
    shells = list()
    for i, (atom, c3d) in enumerate(zip(atoms, coords3d)):
        Zs = str(ATOMIC_NUMBERS[atom.lower()])
        basis_shells = basis[Zs]["electron_shells"]
        for bshell in basis_shells:
            L = bshell["angular_momentum"]
            assert len(L) == 1
            L = L[0]
            exponents = bshell["exponents"]
            for coeffs in bshell["coefficients"]:
                shell = Shell(
                    L=L,
                    center=c3d,
                    coeffs=coeffs,
                    exps=exponents,
                    atomic_num=Zs,
                    center_ind=i,
                )
                shells.append(shell)
    shells = Shells(shells, **kwargs)
    return shells


class Basis:
    """
    Read basis sets from files.
    Bring them in a suitable order. 1s2s2p3s3p4s3d etc.
    """

    pass
