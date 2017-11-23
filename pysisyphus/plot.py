#!/usr/bin/env python3

import argparse
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
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


def load_results(key):
    with open("results.yaml") as handle:
        all_results = yaml.load(handle.read())
    num_cycles = len(all_results)
    print(f"Loaded informations of {num_cycles} cycle(s).")

    tmp_list = list()
    for res_per_cycle in all_results:
            tmp_list.append([res[key] for res in res_per_cycle])
    return np.array(tmp_list), num_cycles


def plot_saras():
    key = "sa_energies"
    sa_ens, num_cycles = load_results(key)
    sa_ens -= sa_ens.min(axis=(2,1), keepdims=True)
    first_cycle = sa_ens[0]

    fig, ax = plt.subplots()
    ax.set_xlabel("Image")
    ax.set_ylabel("ΔE / au")
    lines = ax.plot(first_cycle)

    def animate(i):
        fig.suptitle("Cycle {}".format(i+1))
        for j, line in enumerate(lines):
            line.set_ydata(sa_ens[i][:,j])

    anim = FuncAnimation(
            fig, animate, interval=500, frames=num_cycles
    )

    plt.show()


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
