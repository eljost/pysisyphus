#!/usr/bin/env python3

import argparse
import sys

import matplotlib.pyplot as plt
import pandas as pd

from pysisyphus.constants import AU2KJPERMOL


def plot_energies():
    df = pd.read_csv("energies.csv")
    cmap = plt.get_cmap("Greys")
    df = df.transpose()
    df -= df.values.min()
    df *= AU2KJPERMOL
    ax = df.plot(
            title="Energies",
            colormap=cmap,
    )
    ax.set_xlabel("Image")
    ax.set_ylabel("dE / kJ mol⁻¹")
    plt.show()


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--energies", action="store_true",
                        help="Plot energies.")
    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])
    if args.energies:
        plot_energies()


if __name__ == "__main__":
    run()
