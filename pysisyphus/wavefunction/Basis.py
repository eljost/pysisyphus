import functools
import json
from pathlib import Path
from typing import Tuple

import numpy as np

from pysisyphus.config import BASIS_LIB_DIR
from pysisyphus.elem_data import ATOMIC_NUMBERS
from pysisyphus.Geometry import Geometry
from pysisyphus.wavefunction import Shell, Shells
from pysisyphus.wavefunction.helpers import L_MAP


def basis_from_json(name):
    basis_path = Path(name).with_suffix(".json")
    if not basis_path.is_absolute():
        basis_path = BASIS_LIB_DIR / basis_path

    with open(basis_path) as handle:
        data = json.load(handle)
    elements = data["elements"]
    return elements


def basis_from_orca_str(basis_str):
    """
    %basis
     newgto 1
     s 2
     1 0.7 0.5
     1 0.8 0.6
     p 2
     1 0.7 0.5
     1 0.8 0.6
     d 2
     1 0.7 0.5
     1 0.8 0.6
     d 3
     1 0.7 0.5
     1 0.8 0.6
     1 0.6 0.4
     end
    end
    """
    basis = dict()
    for line in basis_str.strip().split("\n"):
        line = line.strip().lower()
        # Skip comments
        if any([line.startswith(key) for key in ("#", "%basis", "end")]) or line == "":
            continue
        tokens = line.split()
        if tokens[0] == "newgto":
            atomic_num = int(tokens[1])
            basis.setdefault(atomic_num, dict())
        elif tokens[0].isalpha():
            ang_mom, nprims = tokens
            ang_mom = L_MAP[ang_mom]
            nprims = int(nprims)
            prev_counter = 0
            shell = {
                "angular_momentum": ang_mom,
                "coefficients": list(),
                "exponents": list(),
            }
            basis[atomic_num].setdefault("electron_shells", list()).append(shell)
        else:
            counter, exponent, coefficient = [float(t) for t in tokens]
            counter = int(counter)
            if coefficient <= 1e-8:
                continue
            assert counter == (prev_counter + 1)
            shell["coefficients"].append(coefficient)
            shell["exponents"].append(exponent)
    return basis


def basis_from_pyscf_str(basis_str):
    """
    # Comment1
    He    S
         13.6267000              0.1752300
          1.9993500              0.8934830
          0.3829930              0.0000000
    He    S
         13.6267000              0.0000000
          1.9993500              0.0000000
          0.3829930              1.0000000
    """
    basis = dict()
    for line in basis_str.strip().split("\n"):
        line = line.strip().lower()
        # Skip comments
        if line.startswith("#") or (line == "") or ("basis" in line) or ("end" in line):
            continue
        tokens = line.split()
        if tokens[0].isalpha():
            element, ang_mom = [t.lower() for t in tokens]
            ang_mom = L_MAP[ang_mom]
            if ang_mom == "sp":
                raise Exception("SP shells are not yet handled")
            atomic_num = ATOMIC_NUMBERS[element]
            shell = {
                "angular_momentum": ang_mom,
                "coefficients": list(),
                "exponents": list(),
            }
            basis.setdefault(atomic_num, dict()).setdefault(
                "electron_shells", list()
            ).append(shell)
        else:
            exponent, coefficient = [float(t) for t in tokens]
            if coefficient <= 1e-8:
                continue
            shell["coefficients"].append(coefficient)
            shell["exponents"].append(exponent)
    return basis


@functools.singledispatch
def shells_with_basis(
    atoms: Tuple, coords, basis=None, name=None, shells_cls=None, **kwargs
):
    assert (basis is not None) or (name is not None)
    if shells_cls is None:
        shells_cls = Shells
    if name is not None:
        basis = basis_from_json(name)

    coords3d = np.reshape(coords, (len(atoms), 3))
    shells = list()
    for i, (atom, c3d) in enumerate(zip(atoms, coords3d)):
        Zs = str(ATOMIC_NUMBERS[atom.lower()])
        basis_shells = basis[Zs]["electron_shells"]
        for bshell in basis_shells:
            L = bshell["angular_momentum"]
            assert len(L) == 1  # Disallow SP shells for now.
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
    shells = shells_cls(shells, **kwargs)
    return shells


@shells_with_basis.register
def _(geom: Geometry, **kwargs):
    return shells_with_basis(geom.atoms, geom.coords, **kwargs)


class Basis:
    """
    Read basis sets from files.
    Bring them in a suitable order. 1s2s2p3s3p4s3d etc.

    NOT YET USED
    """
    pass
