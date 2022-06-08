import json

import numpy as np

from pysisyphus.helpers_pure import file_or_str
from pysisyphus.wavefunction import get_l, Shell, ORCAShells, Wavefunction
from pysisyphus.wavefunction.helpers import BFType


@file_or_str(".json")
def shells_from_json(text):
    data = json.loads(text)
    atoms = data["Molecule"]["Atoms"]

    _shells = list()
    for atom in atoms:
        center = np.array(atom["Coords"])
        center_ind = atom["Idx"]
        bfs = atom["BasisFunctions"]
        atomic_num = atom["ElementNumber"]
        for bf in bfs:
            exps = np.array(bf["Exponents"])
            coeffs = np.array(bf["Coefficients"])
            L = get_l(bf["Shell"])
            shell = Shell(
                L=L,
                atomic_num=atomic_num,
                center=center,
                coeffs=coeffs,
                exps=exps,
                center_ind=center_ind,
            )
            _shells.append(shell)
    shells = ORCAShells(_shells)
    return shells


@file_or_str(".json")
def wavefunction_from_json(text):
    data = json.loads(text)
    mol = data["Molecule"]

    charge = mol["Charge"]
    mult = mol["Multiplicity"]
    atoms = list()
    coords = list()
    for atom in mol["Atoms"]:
        atom_label = atom["ElementLabel"]
        atom_coords = atom["Coords"]
        atoms.append(atom_label)
        coords.append(atom_coords)

    unrestricted = (mol["HFTyp"] == "UHF") or (mult != 1)
    occ_a = 0
    occ_b = 0
    mos = mol["MolecularOrbitals"]["MOs"]
    Ca = list()
    Cb = list()

    passed_a = False
    prev_occ = -1
    for mo in mos:
        occ = mo["Occupancy"]
        if (occ % 1) != 0.0:
            raise Exception("Fractional occupations are not handled!")
        occ = int(occ)

        if not passed_a and (prev_occ == 0) and (occ == 1):
            passed_a = True

        mo_coeffs = mo["MOCoefficients"]
        if passed_a:
            occ_b += occ
            Cb.append(mo_coeffs)
        else:
            occ_a += occ
            Ca.append(mo_coeffs)
        prev_occ = occ

    # MOs must be in columns
    Ca = np.array(Ca).T
    Cb = np.array(Cb).T
    # Restricted calculation
    if Cb.size == 0:
        Cb = Ca.copy()
        occ_a = occ_b = occ_a // 2
    C = np.stack((Ca, Cb))

    shells = shells_from_json(text)

    return Wavefunction(
        atoms=atoms,
        coords=coords,
        charge=charge,
        mult=mult,
        unrestricted=unrestricted,
        occ=(occ_a, occ_b),
        C=C,
        bf_type=BFType.PURE_SPHERICAL,
        shells=shells,
    )
