# [1] 10.1016/j.comptc.2014.10.002
#     JANPA: An open source cross-platform implementation of the
#     Natural Population Analysis on the Java platform
#     Nikolaienko, Bulavin, Hovorun, 2014

import json

import numpy as np

from pysisyphus.constants import BOHR2ANG
from pysisyphus.helpers_pure import file_or_str
from pysisyphus.io import bson
from pysisyphus.io import molden
from pysisyphus.wavefunction import (
    get_l,
    Shell,
    ORCAShells,
    Wavefunction,
)
from pysisyphus.wavefunction.helpers import BFType


def get_coord_factor(data):
    """Somewhere between ORCA 5.0 and 5.0.4 ORCA switched from Bohr to Angstrom."""
    mol = data["Molecule"]
    unit = mol["CoordinateUnits"]
    if unit == "Angs":
        factor = 1 / BOHR2ANG
    elif unit == "Bohrs":
        factor = 1.0
    else:
        raise Exception(f"Unknown unit '{unit}'!")
    return factor


def shells_from_json_dict(data, **kwargs):
    atoms = data["Molecule"]["Atoms"]
    coord_factor = get_coord_factor(data)

    _shells = list()
    for atom in atoms:
        center = np.array(atom["Coords"]) * coord_factor
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
    shells = ORCAShells(_shells, **kwargs)
    return shells


@file_or_str(".json")
def shells_from_json(text, **kwargs):
    data = json.loads(text)
    return shells_from_json_dict(data, **kwargs)


def wavefunction_from_json_dict(data, **kwargs):
    kwargs = kwargs.copy()
    shell_kwargs = kwargs.pop("shell_kwargs", {})

    mol = data["Molecule"]
    coord_factor = get_coord_factor(data)

    charge = mol["Charge"]
    mult = mol["Multiplicity"]
    atoms = list()
    coords = list()
    for atom in mol["Atoms"]:
        atom_label = atom["ElementLabel"]
        atom_coords = atom["Coords"]
        atoms.append(atom_label)
        coords.append(atom_coords)

    coords = np.array(coords) * coord_factor

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
        mo_ens = list()
        occ = 0
        for mo in mos:
            _occ = mo["Occupancy"]
            if (_occ % 1) != 0.0:
                raise Exception("Fractional occupations are not handled!")
            occ += int(_occ)
            mo_coeffs.append(mo["MOCoefficients"])
            mo_ens.append(mo["OrbitalEnergy"])
        # MOs must be in columns
        mo_coeffs = np.array(mo_coeffs).T
        return occ, mo_coeffs, mo_ens

    occ_a, Ca, ens_a = get_occ_and_mo_coeffs(mos_a)
    occ_b, Cb, ens_b = get_occ_and_mo_coeffs(mos_b)

    # Restricted calculation
    if Cb.size == 0:
        Cb = Ca.copy()
        occ_a = occ_b = occ_a // 2
        ens_b = ens_a
    C = np.stack((Ca, Cb))
    mo_ens = np.stack((ens_a, ens_b))

    shells = shells_from_json_dict(data, **shell_kwargs)

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
        mo_ens=mo_ens,
        **kwargs,
    )


@file_or_str(".json")
def wavefunction_from_json(text, **kwargs):
    data = json.loads(text)
    return wavefunction_from_json_dict(data, **kwargs)


@file_or_str(".bson")
def wavefunction_from_bson(text, **kwargs):
    data = bson.loads(text)
    return wavefunction_from_json_dict(data, **kwargs)


@file_or_str(".molden", ".input")
def wavefunction_from_orca_molden(text, charge=None, **kwargs):
    return molden.wavefunction_from_molden(
        text, charge=charge, orca_contr_renorm=True, **kwargs
    )
