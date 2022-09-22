#!/usr/bin/env python

import argparse
import sys

import numpy as np

from pysisyphus.constants import AU2EV
from pysisyphus.wavefunction import foster_boys


def parse_args(args):
    parser = argparse.ArgumentParser("Diabatization using Foster-Boys-localization.")

    parser.add_argument(
        "ens", type=str, help="File containing adiabatic energies in au."
    )
    parser.add_argument(
        "dpms",
        type=str,
        help="File containing adiabatic dipole moments in au.",
    )
    parser.add_argument(
        "--debye",
        action="store_true",
        help="Adiabatic dipole moments are expected to be given in Debye instead of au.",
    )
    parser.add_argument(
        "--eV",
        action="store_true",
        help="Adiabatc energies are expected to be given in electronvolts instead of Hartree.",
    )
    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    ens_org = np.loadtxt(args.ens, dtype=float)

    assert ens_org.ndim == 1
    nstates = ens_org.size
    
    # Convert energies to eV
    if args.eV:
        en_factor = 1.
        en_unit_org = "eV"
    else:
        en_factor = AU2EV
        en_unit_org = "Eh"
    ens = ens_org.copy()
    ens -= ens[0]
    ens *= en_factor

    print(f"Loaded adiabatic energies of {nstates} states:")
    for i, (en_org, en) in enumerate(zip(ens_org, ens)):
        print(f"{i:02d}: E={en_org: >18.6f} {en_unit_org}\tΔE={en: >18.6f} eV")
    print()

    en_mat = np.diag(ens)
    print("Adiabatic energy Matrix <I|H|J> in eV")
    print(en_mat)
    print()

    dpm_mat = np.zeros((3, nstates, nstates))
    nunique_dpms = sum(range(nstates+1))
    unique_dpms = list()
    if args.debye:
        dpm_fact = 0.393456
        # dpm_unit_org = "Debye"
    else:
        dpm_fact = 1.0
        # dpm_unit_org = "au"

    with open(args.dpms) as handle:
        dpm_text = handle.read()

    for i, line in enumerate(dpm_text.strip().split("\n")):
        floats = [float(_) for _ in line.strip().split()]
        assert len(floats) == 5
        from_, to_, dpmx, dpmy, dpmz = floats
        assert from_.is_integer()
        assert to_.is_integer()
        from_ = int(from_)
        to_ = int(to_)
        dpm_mat[0, from_, to_] = dpmx * dpm_fact
        dpm_mat[1, from_, to_] = dpmy * dpm_fact
        dpm_mat[2, from_, to_] = dpmz * dpm_fact
        unique_dpms.append(frozenset((from_, to_)))
    unique_dpms = set(unique_dpms)
    act_nunique_dpms = len(unique_dpms)
    assert act_nunique_dpms == nunique_dpms
    print(f"Expected {nunique_dpms} unique dipole moments, found {act_nunique_dpms}")

    print("Dipole moment matrices in au")
    print(dpm_mat)
    print()

    U0 = np.eye(nstates)
    dia_result = foster_boys(U0, dpm_mat)
    assert dia_result.is_converged
    U = dia_result.C
    print()
    print("Final rotation matrix U")
    print(U)
    print()

    en_mat_rot = U @ en_mat @ U.T
    ens_dia = np.diag(en_mat_rot)
    print("Diabatic energy matrix (U @ <I|H|J> @ U.T) in eV")
    print(en_mat_rot)
    print()

    print(f"Energies of diabatic states")
    for i, (en) in enumerate(ens_dia):
        print(f"{i:02d}: E={en: >18.6f} eV")
    print()

    print("Diabatic couplings in eV")
    tril = np.triu_indices(nstates, k=1)
    for from_, to_ in zip(*tril):
        coupling = abs(en_mat_rot[from_, to_])
        print(f"{from_:02d} {to_:02d} {coupling: >16.8f} eV")


if __name__ == "__main__":
    run()
