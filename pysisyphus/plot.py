#!/usr/bin/env python3

import argparse
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import splrep, splev
import yaml

from pysisyphus.constants import AU2KJPERMOL
from pysisyphus.peakdetect import peakdetect


def plot_energies(df):
    cmap = plt.get_cmap("Greys")
    df = df.transpose()
    df -= df.values.min()
    df *= AU2KJPERMOL

    fig, ax = plt.subplots()
    ax = df.plot(
            ax=ax,
            title="Energies",
            colormap=cmap,
            legend=True,
    )
    kwargs = {
        "ls": ":",
        "color": "darkgrey",
    }
    try:
        last_row = df.transpose().iloc[-1]
        spl = splrep(last_row.index, last_row)
        images = len(last_row.index)
        # Calculate interpolated values
        x2 = np.linspace(0, images, 100)
        y2 = splev(x2, spl)
        # Only consider maxima
        peak_inds, _ = peakdetect(y2, lookahead=2)
        peak_inds = np.array(peak_inds)[:, 0].astype(int)
        peak_xs = x2[peak_inds]
        peak_ys = y2[peak_inds]
        ax.plot(x2, y2, peak_xs, peak_ys, "x")
        for px, py in zip(peak_xs, peak_ys):
            ax.axhline(y=py, **kwargs)
            line = matplotlib.lines.Line2D([px, px], [0, py], **kwargs)
            ax.add_line(line)
    except TypeError:
        print("Not enough images for splining!")

    # Always draw a line at the minimum y=0
    ax.axhline(y=0, **kwargs)
    ax.set_xlabel("Image")
    ax.set_ylabel("dE / kJ mol⁻¹")
    plt.show()


def load_results(keys=()):
    with open("results.yaml") as handle:
        all_results = yaml.load(handle.read())
    print(f"Loaded informations of {len(all_results)} cycles.")
    if len(keys) == 0:
        return all_results

    tmp_dict = {k: [] for k in keys}
    for res_per_cycle in all_results:
        for k in keys:
            tmp_dict[k].append([res[k] for res in res_per_cycle])

    print(tmp_dict)


def plot_saras():
    keys = ("energy", )
    ens = load_results(keys)


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--energies", action="store_true",
                        help="Plot energies.")
    parser.add_argument("--until", type=int,
                        help="Only show until cycle [until].")
    parser.add_argument("--saras", action="store_true")
    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    if args.energies:
        df = pd.read_csv("energies.csv")
        if args.until:
            df = df[:args.until]
        plot_energies(df)
    elif args.saras:
        plot_saras()


if __name__ == "__main__":
    run()
