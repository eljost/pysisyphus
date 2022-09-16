# [1] 10.1016/j.comptc.2014.10.002
#     JANPA: An open source cross-platform implementation of the
#     Natural Population Analysis on the Java platform
#     Nikolaienko, Bulavin, Hovorun, 2014

import json

import numpy as np
from scipy.special import gamma

from pysisyphus.elem_data import ATOMIC_NUMBERS
from pysisyphus.helpers_pure import file_or_str
from pysisyphus.io.molden import parse_molden
from pysisyphus.wavefunction import get_l, MoldenShells, Shell, ORCAShells, Wavefunction
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
                center=center,
                coeffs=coeffs,
                exps=exps,
                center_ind=center_ind,
                atomic_num=atomic_num,
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
    mos = mol["MolecularOrbitals"]["MOs"]
    if unrestricted:
        mo_num = len(mos)
        assert mo_num % 2 == 0
        a_num = mo_num // 2
        mos_a = mos[:a_num]
        mos_b = mos[a_num:]
    else:  # restricted
        mos_a = mos
        mos_b = list()

    def get_occ_and_mo_coeffs(mos):
        mo_coeffs = list()
        occ = 0
        for mo in mos:
            _occ = mo["Occupancy"]
            if (_occ % 1) != 0.0:
                raise Exception("Fractional occupations are not handled!")
            occ += int(_occ)

            mo_coeffs.append(mo["MOCoefficients"])
        # MOs must be in columns
        mo_coeffs = np.array(mo_coeffs).T
        return occ, mo_coeffs

    occ_a, Ca = get_occ_and_mo_coeffs(mos_a)
    occ_b, Cb = get_occ_and_mo_coeffs(mos_b)

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


def radial_integral(l, exponent):
    """
    Integrates
        (r r**l * exp(-exponent * r**2))**2 dr from r=0 to r=oo
    as described in the SI of the JANPA paper [1] (see top of page 8,
    second integral in the square root.

    In my opinion, the integrals lacks a factor 'r'. Below, some sympy code
    can be found to solve this integral (including 1*r).

    import sympy as sym
    r, z = sym.symbols("r z", positive=True)
    l = sym.symbols("l", integer=True, positive=True)
    sym.integrate((r * r**l * sym.exp(-z*r**2))**2, (r, 0, sym.oo))

    The 'solved' integral on page 8 is correct again.

     ∞
     ⌠
     ⎮              2
     ⎮  2  2⋅l  -2⋅r ⋅z
     ⎮ r ⋅r   ⋅ℯ        dr = (2*z)**(-l - 1/2)*gamma(l + 3/2)/(4*z)
     ⌡
     0
    """
    return (2 * exponent) ** (-l - 1 / 2) * gamma(l + 3 / 2) / (4 * exponent)


@file_or_str(".molden", ".input")
def shells_from_molden(text):
    data = parse_molden(text, with_mos=False)

    atoms, coords, _ = parse_molden_atoms(data)
    atomic_numbers = [ATOMIC_NUMBERS[atom.lower()] for atom in atoms]
    coords3d = coords.reshape(-1, 3)

    dividers = {
        0: 1,
        1: 1,
        2: 3,
        3: 15,
        4: 35,
    }

    def fix_contr_coeffs(l, coeffs, exponents):
        """Fix contraction coefficients. Based on equations found in the SI
        of the JANPA paper [1]."""
        l = get_l(l)
        divider = dividers[l]
        rad_ints = radial_integral(l, exponents)
        norms2 = rad_ints * 4 * np.pi / (2 * l + 1) / divider
        norms = np.sqrt(norms2)
        normed_coeffs = coeffs * norms
        return normed_coeffs

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

    shells = MoldenShells(_shells)
    return shells


@file_or_str(".molden", ".input")
def wavefunction_from_molden(text, charge=None, **wf_kwargs):
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

    molden_charge = nuc_charge - (occ_a + occ_b)
    if charge is None:
        charge = molden_charge
    # Multiplicity = (2S + 1), S = number of unpaired elecs. * 0.5
    mult = (occ_a - occ_b) + 1
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
        **wf_kwargs,
    )
