import argparse
import sys

from pysisyphus.config import p_DEFAULT, T_DEFAULT
from pysisyphus.thermo import get_thermoanalysis_from_hess_h5, print_thermoanalysis


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("hess_h5", help="HDF5 Hessian file from pysisyphus.")
    parser.add_argument("-T", default=T_DEFAULT, type=float, help="Temperature")
    parser.add_argument("-p", default=p_DEFAULT, type=float, help="Pressure")
    parser.add_argument("--pg", default="c1", help="Point group.")
    parser.add_argument(
        "--calorie",
        action="store_true",
        help="Output in kcal mol⁻¹ instead of kJ mol⁻¹.",
    )

    return parser.parse_args(args)


def run_thermo():
    args = parse_args(sys.argv[1:])
    hess_h5 = args.hess_h5
    T = args.T
    p = args.p
    point_group = args.pg
    unit = "calorie" if args.calorie else "joule"
    thermo = get_thermoanalysis_from_hess_h5(hess_h5, T=T, p=p, point_group=point_group)
    print_thermoanalysis(thermo, unit=unit)
