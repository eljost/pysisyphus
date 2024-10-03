#!/usr/bin/env python3

import argparse
import sys

import numpy as np

from pysisyphus.Geometry import Geometry
from pysisyphus.io import hessian


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("atom", type=str)
    parser.add_argument("energy", type=float)
    parser.add_argument("mult", type=int)
    parser.add_argument("--charge", default=0, type=int)
    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    atom = args.atom
    energy = args.energy
    mult = args.mult
    charge = args.charge
    out_fn = f"{atom}_final_hessian.h5"

    geom = Geometry(atom, np.zeros(3))
    hessian.save_hessian(
        out_fn,
        geom,
        cart_hessian=np.zeros((3, 3)),
        energy=energy,
        mult=mult,
        charge=charge,
    )
    print(f"Saved Hessian to '{out_fn}'.")


if __name__ == "__main__":
    run()
