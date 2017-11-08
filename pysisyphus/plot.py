#!/usr/bin/env python3

import argparse
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import splrep, splev

from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.peakdetect import peakdetect


def plot_energies(df):
    cmap = plt.get_cmap("Greys")
    df = df.transpose()
    df -= df.values.min()
    df *= AU2KJPERMOL

    last_row = df.transpose().iloc[-1]
    spl = splrep(last_row.index, last_row)
    images = len(last_row.index)

    # Calculate interpolated values
    x2 = np.linspace(0, images, 100)
    y2 = splev(x2, spl)
    # Only consider maxima
    peak_inds, _ = peakdetect(y2, lookahead=2)
    peak_inds = np.array(peak_inds)[:,0].astype(int)
    peak_xs = x2[peak_inds]
    peak_ys = y2[peak_inds]

    fig, ax = plt.subplots()
    ax = df.plot(
            ax=ax,
            title="Energies",
            colormap=cmap,
            legend=False,
    )
    ax.plot(x2, y2, peak_xs, peak_ys, "x")
    kwargs = {
        "ls": ":",
        "color": "darkgrey",
    }
    # Always draw a line at the minimum y=0
    ax.axhline(y=0, **kwargs)
    for px, py in zip(peak_xs, peak_ys):
        ax.axhline(y=py, **kwargs)
        line = matplotlib.lines.Line2D([px, px], [0, py], **kwargs)
        ax.add_line(line)
    ax.set_xlabel("Image")
    ax.set_ylabel("dE / kJ mol⁻¹")
    plt.show()


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--energies", action="store_true",
                        help="Plot energies.")
    parser.add_argument("--until", type=int, default=-1,
                        help="Only show until cycle [until].")
    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])
    if args.energies:
        df = pd.read_csv("energies.csv")
        plot_energies(df.head(args.until))


if __name__ == "__main__":
    run()
