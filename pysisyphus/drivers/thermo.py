import argparse
import sys

from pysisyphus.thermo import get_thermoanalysis_from_hess_h5, print_thermoanalysis


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("hess_h5", help="HDF5 Hessian file from pysisyphus.")
    parser.add_argument("--T", default=298.15, type=float, help="Temperature")
    parser.add_argument("--pg", default="c1", help="Point group.")

    return parser.parse_args(args)


def run_thermo():
    args = parse_args(sys.argv[1:])
    hess_h5 = args.hess_h5
    T = args.T
    point_group = args.pg
    thermo = get_thermoanalysis_from_hess_h5(hess_h5, T=T, point_group=point_group)
    print_thermoanalysis(thermo)
