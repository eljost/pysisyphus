#!/usr/bin/env python3

import argparse
from pathlib import Path
import sys

import numpy as np

from pysisyphus.calculators.ORCA import parse_orca_cis, parse_orca_all_energies
from pysisyphus.calculators.Gaussian16 import EXC_STATE_RE, parse_ci_coeffs
from pysisyphus.constants import AU2EV
from pysisyphus.wavefunction import Wavefunction
from pysisyphus.wavefunction import so_coupling


def run_orca(base: Path):

    wf_fn = base.with_suffix(".bson")
    cis_fn = base.with_suffix(".cis")
    log_fn = base.with_suffix(".log")

    wf = Wavefunction.from_file(wf_fn)
    Xa, Ya, Xb, Yb = parse_orca_cis(
        cis_fn, restricted_same_ab=True, triplets_only=False
    )

    sing_ens = parse_orca_all_energies(log_fn, triplets=False, do_tddft=True)
    # Drop GS from triplets
    trip_ens = parse_orca_all_energies(log_fn, triplets=True, do_tddft=True)[1:]
    nroots = len(sing_ens) - 1
    singlets = list(range(nroots))
    triplets = list(range(nroots, 2 * nroots))
    Xas = Xa[singlets]
    Yas = Ya[singlets]
    Xat = Xa[triplets]
    Yat = Ya[triplets]

    gs_en = sing_ens[0]
    sing_ens = sing_ens - gs_en
    trip_ens = trip_ens - gs_en

    socs = so_coupling.run(wf, Xas, Yas, Xat, Yat, sing_ens, trip_ens)
    return socs


def run_gaussian(base: Path):

    wf_fn = base.with_suffix(".fchk")
    log_fn = base.with_suffix(".log")

    wf = Wavefunction.from_file(wf_fn)
    Xa, Ya, Xb, Yb = parse_ci_coeffs(log_fn)

    with open(log_fn) as handle:
        text = handle.read()
    exc_states = EXC_STATE_RE.findall(text)

    singlets = list()
    triplets = list()
    sing_ens = [0.0]
    trip_ens = []
    for state in exc_states:
        ind, spin, exc_en_eV, *_ = state
        ind = int(ind) - 1
        if "Triplet" in spin:
            triplets.append(ind)
            trip_ens.append(exc_en_eV)
        elif "Singlet" in spin:
            singlets.append(ind)
            sing_ens.append(exc_en_eV)
        else:
            raise Exception("Unknown multiplicity!")
    sing_ens = np.array(sing_ens, dtype=float) / AU2EV
    trip_ens = np.array(trip_ens, dtype=float) / AU2EV

    Xas = Xa[singlets]
    Yas = Ya[singlets]
    Xat = Xa[triplets]
    Yat = Ya[triplets]

    socs = so_coupling.run(wf, Xas, Yas, Xat, Yat, sing_ens, trip_ens)
    return socs


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("base")
    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    base = Path(args.base)

    log = base.with_suffix(".log")
    with open(log) as handle:
        text = handle.read()

    if "O   R   C   A" in text:
        run_orca(base)
    elif "Entering Gaussian System" in text:
        run_gaussian(base)
    else:
        raise Exception("Unknown program!")


if __name__ == "__main__":
    run()
