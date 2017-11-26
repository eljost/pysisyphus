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

from pysisyphus.constants import AU2KJPERMOL, BOHR2ANG
from pysisyphus.peakdetect import peakdetect



class Plotter:
    def __init__(self, data, ylabel, interval=750):
        self.data = data
        self.ylabel = ylabel
        self.interval = interval

        self.cycles = len(self.data)
        self.pause = True

        self.fig, self.ax = plt.subplots()
        self.fig.canvas.mpl_connect('key_press_event', self.on_keypress)

        self.ax.set_xlabel("Image")
        y_min = self.data.min()#*1.05
        y_max = self.data.max()#*1.05
        self.ax.set_ylim(y_min, y_max)
        self.ax.set_ylabel(self.ylabel)

        if self.data.ndim == 2:
            self.update_func = self.update_plot
        elif self.data.ndim == 3:
            self.update_func = self.update_plot2

    def update_plot(self, i):
        self.fig.suptitle("Cycle {}".format(i+1))
        self.lines[0].set_ydata(self.data[i])

    def update_plot2(self, i):
        self.fig.suptitle("Cycle {}".format(i+1))
        for j, line in enumerate(self.lines):
            line.set_ydata(self.data[i][:,j])

    def animate(self):
        self.lines = self.ax.plot(self.data[0], "o-")
        self.animation = matplotlib.animation.FuncAnimation(self.fig,
                                                            self.update_func,
                                                            frames=self.cycles,
                                                            interval=self.interval)
        plt.show()

    def on_keypress(self, event):
        """Pause on SPACE press."""
        #https://stackoverflow.com/questions/41557578
        if event.key == " ":
            if self.pause:
                self.animation.event_source.stop()
            else:
                self.animation.event_source.start()
            self.pause = not self.pause


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

    plotter = Plotter(sa_ens, "ΔE / au")
    plotter.animate()


def plot_params(inds):
    def get_bond_length(coords_slice):
        return np.linalg.norm(coords_slice[0]-coords_slice[1]) * BOHR2ANG * 100

    def get_angle(coords_slice):
        vec1 = coords_slice[0] - coords_slice[1]
        vec2 = coords_slice[2] - coords_slice[1]
        vec1n = np.linalg.norm(vec1)
        vec2n = np.linalg.norm(vec2)
        dotp = np.dot(vec1, vec2)
        radians = np.arccos(dotp / (vec1n * vec2n))
        return radians * 180 / np.pi

    type_dict = {
        2: ("bond length / pm", get_bond_length),
        3: ("angle / °", get_angle),
    }
    ylabel, func = type_dict[len(inds)]

    key = "coords"
    coords, num_cycles = load_results(key)
    ac = list()
    for i, per_cycle in enumerate(coords):
        pc = list()
        for j, per_image in enumerate(per_cycle):
            coords_slice = per_image.reshape(-1, 3)[inds]
            param = func(coords_slice)
            pc.append(param)
        ac.append(pc)
    ac_arr = np.array(ac)

    plotter = Plotter(ac_arr, ylabel)
    plotter.animate()


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--energies", action="store_true",
                        help="Plot energies.")
    parser.add_argument("--until", type=int,
                        help="Only show until cycle [until].")
    parser.add_argument("--saras", action="store_true")
    parser.add_argument("--params", type=int, nargs="+")
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
    elif args.params:
        plot_params(args.params)


if __name__ == "__main__":
    run()
