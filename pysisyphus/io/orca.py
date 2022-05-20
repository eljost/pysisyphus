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

    unrestricted = mol["HFTyp"] == "UHF"
    occ = 0
    mos = mol["MolecularOrbitals"]["MOs"]
    C = list()
    for mo in mos:
        occ += mo["Occupancy"]
        C.append(mo["MOCoefficients"])
    C = np.array(C)
    if unrestricted:
        C = C.reshape(2, C.shape[-1], -1)
        C = np.transpose(C, (0, 2, 1))
    else:
        C = C.T
    occ_a = occ_b = int(occ) // 2

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
