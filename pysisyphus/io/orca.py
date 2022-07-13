import json

import numpy as np
from scipy.special import factorial2

from pysisyphus.elem_data import ATOMIC_NUMBERS
from pysisyphus.helpers_pure import file_or_str
from pysisyphus.io.molden import parse_molden
from pysisyphus.wavefunction import get_l, Shell, ORCAShells, Wavefunction
from pysisyphus.wavefunction.helpers import BFType, get_l


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


def parse_molden_atoms(data):
    atoms = list()
    coords = list()
    nuc_charges = list()
    for atom in data["atoms"]:
        atoms.append(atom["symbol"])
        coords.extend(atom["xyz"])
        Z = atom["atomic_number"]
        nuc_charges.append(Z)
    atoms = tuple(atoms)
    coords = np.array(coords, dtype=float).flatten()
    nuc_charges = np.array(nuc_charges, dtype=int)
    return atoms, coords, nuc_charges


@file_or_str(".molden", ".input")
def shells_from_molden(text):
    data = parse_molden(text, with_mos=False)

    atoms, coords, _ = parse_molden_atoms(data)
    atomic_numbers = [ATOMIC_NUMBERS[atom.lower()] for atom in atoms]
    coords3d = coords.reshape(-1, 3)

    def fix_contr_coeffs(l, coeffs, exponents):
        l = get_l(l)
        N = (
            (2 ** (l + 0.75) * exponents ** (l / 2 + 0.75))
            / np.pi ** 0.75
            / np.sqrt(factorial2(2 * l - 1))
        )
        fixed_coeffs = coeffs / N
        return fixed_coeffs

    _shells = list()
    for (atom_gtos, center, atomic_num) in zip(data["gto"], coords3d, atomic_numbers):
        center_ind = atom_gtos["number"] - 1
        for shell in atom_gtos["shells"]:
            L = shell["ang_mom"]
            exps = shell["exponents"][0]
            coeffs = shell["coeffs"][0]
            fixed_coeffs = fix_contr_coeffs(L, coeffs, exps)
            shell = Shell(
                L=L,
                center=center,
                coeffs=fixed_coeffs,
                exps=exps,
                center_ind=center_ind,
                atomic_num=atomic_num,
            )
            _shells.append(shell)
    from pysisyphus.wavefunction.shells import ORCAMoldenShells
    # shells = ORCAShells(_shells)
    shells = ORCAMoldenShells(_shells)
    return shells


@file_or_str(".molden", ".input")
def wavefunction_from_molden(text):
    data = parse_molden(text)
    atoms, coords, nuc_charges = parse_molden_atoms(data)
    nuc_charge = sum(nuc_charges)

    spins = list()
    occ_a = 0.0
    occ_b = 0.0
    Ca = list()
    Cb = list()
    for mo in data["mos"]:
        spin = mo["spin"]
        occ = mo["occup"]
        spins.append(spin)
        coeffs = mo["coeffs"]
        if spin == "Alpha":
            occ_a += occ
            Ca.append(coeffs)
        elif spin == "Beta":
            occ_b += occ
            Cb.append(coeffs)
        else:
            raise Exception(f"Spin can only be 'Alpha' or 'Beta', but got '{spin}'!")
    assert occ_a.is_integer
    assert occ_b.is_integer
    occ_a = int(occ_a)
    occ_b = int(occ_b)

    # MOs must be in columns
    Ca = np.array(Ca).T
    Cb = np.array(Cb).T
    # Restricted calculation
    if Cb.size == 0:
        Cb = Ca.copy()
        occ_a = occ_b = occ_a // 2
    C = np.stack((Ca, Cb))

    charge = nuc_charge - (occ_a + occ_b)
    mult = (occ_a - occ_b) * 2 + 1
    unrestricted = set(spins) == {"Alpha", "Beta"}

    shells = shells_from_molden(text)

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
