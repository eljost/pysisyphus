#!/usr/bin/env python3

import argparse
import re
import sys

import numpy as np

from pysisyphus.calculators import ORCA
from pysisyphus.calculators.ORCA import make_sym_mat
from pysisyphus.Geometry import Geometry
from pysisyphus.io.hessian import save_hessian


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("log", help="ORCA log file.")
    parser.add_argument("hess", help="ORCA .hess file.")
    parser.add_argument(
        "--out", default="orca_hessian.h5", help="Pysisyphus HDF5 Hessian output file."
    )

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    log_fn = args.log
    hess_fn = args.hess
    h5_fn = args.out

    parsed = ORCA.parse_hess_file(hess_fn)
    print(f"Read '{hess_fn}'.")
    cart_hessian = make_sym_mat(parsed["hessian"])
    atoms, _, coords3d = zip(*parsed["atoms"][2:])
    coords3d = np.array([c.asList() for c in coords3d])
    geom = Geometry(atoms, coords3d)

    with open(log_fn) as handle:
        log_text = handle.read()
    print(f"Read '{log_fn}'")
    energies = re.findall(r"FINAL SINGLE POINT ENERGY\s+([\d\-\.]+)", log_text)
    energy = float(energies[-1])
    print(
        f"Found {len(energies)} energies in '{log_fn}'.\n"
        f"Using last one: {energy:.6f} au."
    )
    mult_re = re.compile("Multiplicity\s+Mult\s+\.{4}\s+(\d+)")
    mult = int(mult_re.search(log_text).group(1))
    print(f"Multiplicity: {mult}")

    save_hessian(h5_fn, geom, cart_hessian=cart_hessian, energy=energy, mult=mult)
    print(f"Wrote pysisyphus HDF5 Hessian to '{h5_fn}'")


if __name__ == "__main__":
    run()
