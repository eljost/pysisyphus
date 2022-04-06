#!/usr/bin/env python3

import argparse
from pathlib import Path
import shutil
import subprocess
import struct
import sys

import numpy as np

from pysisyphus.wrapper.mwfn import get_mwfn_exc_str, make_cdd


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="Wrapper script for "
        "charge density difference (CDD) cube generation using "
        "(ORCA), pysisyphus and Multiwfn. The orca_2mkl tool is required "
        "when this script is called with a .gbw file. Pysisyphus and Multiwfn "
        "are required always. Currently, only the X-vector is used for CDD generation."
    )
    parser.add_argument("wfn", help="ORCA .gbw or molden file.")
    parser.add_argument("cis", help="ORCA .cis file.")
    parser.add_argument(
        "states", type=int, nargs="+", help="State indices, for which to create CDDs."
    )
    return parser.parse_args(args)


def gbw2molden(gbw):
    pgbw = Path(gbw)
    no_ext = pgbw.stem
    cmd = f"orca_2mkl {no_ext} orca -molden".split()
    _ = subprocess.check_output(
        cmd,
        universal_newlines=True,
    )
    out_fn = pgbw.with_suffix(".molden.input")
    molden_fn = pgbw.with_suffix(".molden")
    shutil.move(out_fn, molden_fn)
    return molden_fn


def parse_cis(cis):
    """
    Read binary CI vector file from ORCA.
        Adapted from TheoDORE 1.7.1, Authors: S. Mai, F. Plasser
        https://sourceforge.net/p/theodore-qc
    """
    energies = list()

    cis_handle = open(cis, "rb")
    # self.log(f"Parsing CI vectors from {cis_handle}")

    # the header consists of 9 4-byte integers, the first 5
    # of which give useful info.
    nvec = struct.unpack("i", cis_handle.read(4))[0]
    # header array contains:
    # [0] index of first alpha occ,  is equal to number of frozen alphas
    # [1] index of last  alpha occ
    # [2] index of first alpha virt
    # [3] index of last  alpha virt, header[3]+1 is equal to number of bfs
    # [4] index of first beta  occ,  for restricted equal to -1
    # [5] index of last  beta  occ,  for restricted equal to -1
    # [6] index of first beta  virt, for restricted equal to -1
    # [7] index of last  beta  virt, for restricted equal to -1
    header = [struct.unpack("i", cis_handle.read(4))[0] for i in range(8)]

    # Assert that all flags regarding unrestricted calculations are -1
    if any([flag != -1 for flag in header[4:8]]):
        raise Exception("parse_cis, no support for unrestricted MOs")

    nfrzc = header[0]
    nocc = header[1] + 1
    nact = nocc - nfrzc
    nmo = header[3] + 1
    nvir = nmo - header[2]
    lenci = nact * nvir
    # self.log(f"nmo = {nmo}, nocc = {nocc}, nact = {nact}, nvir = {nvir}")

    # Loop over states. For non-TDA order is: X+Y of 1, X-Y of 1,
    # X+Y of 2, X-Y of 2, ...
    prev_root = -1
    prev_mult = 1
    iroot_triplets = 0

    # Flags that may later be set to True
    triplets = False
    tda = False
    Xs = list()
    Ys = list()

    for ivec in range(nvec):
        # header of each vector
        # contains 6 4-byte ints, then 1 8-byte double, then 8 byte unknown
        nele, d1, mult, d2, iroot, d3 = struct.unpack("iiiiii", cis_handle.read(24))

        # Will evaluate True only once when triplets were requested.
        if prev_mult != mult:
            triplets = True
            prev_root = -1

        # When we encounter the second "state" we can decide if it is a TDA
        # calculation (without Y-vector).
        if (ivec == 1) and (iroot == prev_root + 1):
            tda = True

        if triplets:
            iroot = iroot_triplets

        ene, d3 = struct.unpack("dd", cis_handle.read(16))
        energies.append(ene)
        # self.log(f"ivec={ivec}, nele={nele}, mult={mult}, iroot={iroot}")
        # Then come nact * nvirt 8-byte doubles with the coefficients
        coeffs = struct.unpack(lenci * "d", cis_handle.read(lenci * 8))
        coeffs = np.array(coeffs).reshape(-1, nvir)
        # create full array, i.e nocc x nvirt
        coeffs_full = np.zeros((nocc, nvir))
        coeffs_full[nfrzc:] = coeffs

        # In this case, we have a non-TDA state, where Y is present!
        # We can recover the original X and Y by first computing X as
        #   X = (X+Y + X-Y) / 2
        # and then
        #   Y = X+Y - X
        if prev_root == iroot:
            X_plus_Y = Xs[-1]
            X_minus_Y = coeffs_full
            X = 0.5 * (X_plus_Y + X_minus_Y)
            Y = X_plus_Y - X
            Xs[-1] = X
            Ys[-1] = Y
        else:
            Xs.append(coeffs_full)
            Ys.append(np.zeros_like(coeffs_full))

        # Somehow ORCA stops to update iroot correctly after the singlet states.
        if (mult == 3) and (tda or (ivec % 2) == 1):
            iroot_triplets += 1

        prev_root = iroot
        prev_mult = mult
    cis_handle.close()

    Xs = np.array(Xs)
    Ys = np.array(Ys)

    # Only return triplet states if present
    if triplets:
        assert (len(Xs) % 2) == 0
        states = len(Xs) // 2
        Xs = Xs[states:]
        Ys = Ys[states:]
    return Xs, Ys, energies


def run():
    args = parse_args(sys.argv[1:])

    wfn = args.wfn
    cis = args.cis
    states = args.states

    if wfn.endswith(".gbw"):
        wfn = gbw2molden(wfn)

    Xs, Ys, energies = parse_cis(cis)
    energies = np.array(energies)

    exc_str = get_mwfn_exc_str(energies, Xs)
    exc_fn = "exc"
    with open(exc_fn, "w") as handle:
        handle.write(exc_str)

    cubes = list()
    for state in states:
        cube = make_cdd(wfn, state, exc_fn)[0]
        cubes.append(cube)
        print(f"Created '{cube}'.")


if __name__ == "__main__":
    run()
