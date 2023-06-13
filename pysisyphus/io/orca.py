# [1] 10.1016/j.comptc.2014.10.002
#     JANPA: An open source cross-platform implementation of the
#     Natural Population Analysis on the Java platform
#     Nikolaienko, Bulavin, Hovorun, 2014

import json

import numpy as np
from scipy.special import gamma

from pysisyphus.helpers_pure import file_or_str
from pysisyphus.io import molden
from pysisyphus.wavefunction import (
    get_l,
    Shell,
    ORCAShells,
    ORCAMoldenShells,
    Wavefunction,
)
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
    molden_shells = molden.shells_from_molden(text)

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
    for shell in molden_shells.shells:
        L, center, _, exps, *_ = shell.as_tuple()
        fixed_coeffs = fix_contr_coeffs(L, shell.coeffs_org, exps)
        fixed_shell = Shell(
            L=L,
            center=center,
            coeffs=fixed_coeffs,
            exps=exps,
            center_ind=shell.center_ind,
            atomic_num=shell.atomic_num,
        )
        _shells.append(fixed_shell)
    shells = ORCAMoldenShells(_shells)
    return shells


@file_or_str(".molden", ".input")
def wavefunction_from_molden(text, charge=None, **wf_kwargs):
    shells_func = shells_from_molden
    return molden.wavefunction_from_molden(
        text, charge=charge, shells_func=shells_func, **wf_kwargs
    )
