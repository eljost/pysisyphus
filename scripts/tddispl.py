import argparse
import sys

from pysisyphus.irc.initial_displ import cubic_displ_for_h5


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("h5_fn")
    parser.add_argument("--dE", type=float, default=-5e-4)

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    h5_fn = args.h5_fn
    dE = args.dE

    step = cubic_displ_for_h5(h5_fn, dE=dE)
    print(step)


if __name__ == "__main__":
    run()
