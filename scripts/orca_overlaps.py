#!/usr/bin/env python3

import argparse
import sys

from pysisyphus.calculators import ORCA
from pysisyphus.helpers_pure import highlight_text


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("--logs", nargs="+", required=True)
    parser.add_argument("--gbws", nargs="+", required=True)
    parser.add_argument("--ciss", nargs="+", required=True)

    parser.add_argument("--h5", default="overlap_data.h5")

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    print(highlight_text("ORCA OVERLAP DATA DUMPER"), "\n")

    logs = args.logs
    gbws = args.gbws
    ciss = args.ciss

    llogs = len(logs)
    lgbws = len(gbws)
    lciss = len(ciss)

    assert llogs == lgbws == lciss, f"Different number of {llogs=}, {lgbws=}, {lciss=}!"

    calc = ORCA(keywords="", dump_fn=args.h5)
    calc.do_tddft = True

    for i, (log, gbw, cis) in enumerate(zip(logs, gbws, ciss)):
        print(f"Processing step {i:03d}:\n\t{log=}\n\t{gbw=}\n\t{cis=}")
        root, triplets = calc.parse_engrad_info(log)
        if root is None:
            root = 1
            print(f"\tNo iroot found for '{log=}! Setting {root=}.")
        calc.root = root
        atoms, coords = calc.parse_atoms_coords(log)
        mo_coeffs, _ = calc.parse_gbw(gbw)
        X, Y = calc.parse_cis(cis)
        with open(log) as handle:
            text = handle.read()
        all_energies = calc.parse_all_energies(text, triplets=triplets)
        overlap_data = (mo_coeffs, X, Y, all_energies)
        calc.store_overlap_data(atoms, coords, path=None, overlap_data=overlap_data)
        print(f"\tParsed step {i:03d}. Calculated root is {root}")
        print()
    calc.dump_overlap_data()
    print(f"Dumped HDF5 to {calc.dump_fn=}")
    print()


if __name__ == "__main__":
    run()
