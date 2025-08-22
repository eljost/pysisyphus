#!/usr/bin/env python


import argparse
import itertools as it
from pathlib import Path
import subprocess
import sys
from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

from pysisyphus.calculators.ORCA import parse_orca_all_energies, parse_orca_cis
from pysisyphus.constants import AU2EV
from pysisyphus.intcoords.setup import get_fragments
from pysisyphus.drivers.spectrum_ctnum import ct_number_plot
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction.excited_states import norm_ci_coeffs, ct_numbers_for_states


def load_data(
    base_name: Path, triplets: bool = False, ignore_bonds: Optional[list[int]] = None
):
    if ignore_bonds is None:
        ignore_bonds = list()

    dump_fn = base_name.with_suffix(".npz")

    if dump_fn.exists():
        data = np.load(dump_fn)
        print(f"Loaded previously calculated data from '{dump_fn}'.")
    else:
        cis_fn = base_name.with_suffix(".cis")
        log_fn = base_name.with_suffix(".log")
        wf_fn = base_name.with_suffix(".bson")

        if not log_fn.exists():
            log_fn = base_name.with_suffix(".out")

        if not wf_fn.exists():
            print(f"{wf_fn} does not exist. Trying to create it with orca_2json")
            cmd = f"orca_2json {log_fn.stem} -bson".split()
            subprocess.run(cmd)

        # Drop GS, only keep excitation energies
        all_ens = parse_orca_all_energies(log_fn, triplets=triplets, do_tddft=True)
        exc_ens = all_ens[1:] - all_ens[0]
        print(f"Parsed excitation energies from '{log_fn}'.")

        # TODO: or create this in the script
        # Load wavefunction
        wf = Wavefunction.from_file(wf_fn)
        print(f"Loaded wavefunction from '{wf_fn}'.")

        frags = get_fragments(
            wf.atoms, wf.coords, ignore_bonds=ignore_bonds, with_unconnected_atoms=True
        )
        # Convert from list of frozensets to list of lists
        frags = list(map(list, frags))

        Xa, Ya, Xb, Yb = parse_orca_cis(cis_fn, restricted_same_ab=True)
        Xa, Ya, Xb, Yb = norm_ci_coeffs(Xa, Ya, Xb, Yb)
        print(f"Parsed CI-coefficients energies from '{cis_fn}'.")

        ct_numbers = ct_numbers_for_states(wf, frags, Xa + Ya, Xb + Yb)
        print("Calculated charge-transfer numbers.")

        tdms = wf.get_transition_dipole_moment(Xa + Ya, Xb + Yb)
        print(f"Calculated transition dipole moments.")
        foscs = wf.oscillator_strength(exc_ens, tdms)
        print(f"Calculated oscillator strengths.")

        homogeneous_frags = np.zeros((len(wf.atoms), 2), dtype=int)
        i = 0
        for j, frag in enumerate(frags):
            for atom in frag:
                homogeneous_frags[i] = j, atom
                i += 1
        data = dict(
            atoms=wf.atoms,
            coords=wf.coords,
            all_ens=all_ens,
            exc_ens=exc_ens,
            Xa=Xa,
            Ya=Ya,
            Xb=Xb,
            Yb=Yb,
            tdms=tdms,
            foscs=foscs,
            ct_numbers=ct_numbers,
            frags=homogeneous_frags,
        )
        np.savez(dump_fn, **data)
    return data


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("base_name")
    parser.add_argument("--triplets", action="store_true")
    parser.add_argument("--bonds", nargs="+", type=int)

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    base_name = Path(args.base_name)
    triplets = args.triplets
    bonds = args.bonds
    if bonds is None:
        bonds = list()

    assert len(bonds) % 2 == 0, "Length of bonds must be a multiple of 2!"
    nbonds = len(bonds) // 2
    ignore_bonds = list()
    for i in range(nbonds):
        ignore_bonds.append(bonds[2 * i : 2 * (i + 1)])

    data = load_data(base_name, triplets, ignore_bonds)

    atoms = data["atoms"]
    coords = data["coords"]
    all_energies = data["all_ens"]
    foscs = data["foscs"]
    ct_numbers = data["ct_numbers"]
    homogenous_frags = data["frags"]
    nfrags = len(homogenous_frags)
    exc_ens = all_energies[1:] - all_energies[0]
    exc_ens_eV = exc_ens * AU2EV
    for i, exc_en in enumerate(exc_ens_eV):
        ct_sum = ct_numbers[i].sum()
        print(f"State {i+1:03d}: {exc_en: >8.2f} eV, ΣCTnums = {ct_sum: >8.4f}")

    frags = []
    for key, _frag in it.groupby(
        homogenous_frags, key=lambda frag_atom_ind: frag_atom_ind[0]
    ):
        cur_frag = []
        for _, atom_ind in _frag:
            cur_frag.append(atom_ind)
        frags.append(cur_frag)

    print(f"Found {nfrags} fragments")
    for i, frag in enumerate(frags):
        print(f"Fragment {i:03d}: {len(frag)} atoms")

    fig = ct_number_plot(atoms, coords, all_energies, foscs, ct_numbers, frags=frags)
    plt.show()


if __name__ == "__main__":
    run()
