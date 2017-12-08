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
from pysisyphus.cos.NEB import NEB
from pysisyphus.Geometry import Geometry
from pysisyphus.peakdetect import peakdetect


class Plotter:
    def __init__(self, coords, data, ylabel, interval=750):
        self.coords = coords
        self.data = data
        self.ylabel = ylabel
        self.interval = interval

        self.cycles = len(self.data)
        self.pause = True

        self.fig, self.ax = plt.subplots()
        self.fig.canvas.mpl_connect('key_press_event', self.on_keypress)

        self.ax.set_xlabel("Path length / Bohr")
        y_min = self.data.min()
        y_max = self.data.max()
        self.ax.set_ylim(y_min, y_max)
        self.ax.set_ylabel(self.ylabel)

        if self.data.ndim == 2:
            self.update_func = self.update_plot
        elif self.data.ndim == 3:
            self.update_func = self.update_plot2

        self.coord_diffs = self.get_coord_diffs(self.coords)

    def get_coord_diffs(self, coords, normalize=False):
        coord_diffs = list()
        for per_cycle in coords:
            tmp_list = [0, ]
            for i in range(len(per_cycle)-1):
                diff = np.linalg.norm(per_cycle[i+1]-per_cycle[i])
                tmp_list.append(diff)
            tmp_list = np.cumsum(tmp_list)
            if normalize:
                tmp_list /= tmp_list.max()
            coord_diffs.append(tmp_list)
        return np.array(coord_diffs)

    def update_plot(self, i):
        self.fig.suptitle("Cycle {}".format(i+1))
        self.lines[0].set_ydata(self.data[i])

    def update_plot2(self, i):
        self.fig.suptitle("Cycle {}".format(i+1))
        for j, line in enumerate(self.lines):
            line.set_ydata(self.data[i][:, j])

    def animate(self):
        self.lines = self.ax.plot(self.coord_diffs[0], self.data[0], "o-")
        self.animation = matplotlib.animation.FuncAnimation(
                                                self.fig,
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


def plot_energies():
    keys = ("energy", "coords")
    (energies, coords), num_cycles, num_images = load_results(keys)

    df = pd.DataFrame(energies)
    cmap = plt.get_cmap("Greys")
    df = df.transpose()
    df -= df.values.min()
    df *= AU2KJPERMOL

    fig, ax = plt.subplots()
    ax = df.plot(
            ax=ax,
            title="Energies",
            colormap=cmap,
            legend=False,
            marker="o",
            xticks=range(num_images),
            xlim=(0, num_images-1),
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
        if not peak_inds:
            ax.plot(x2, y2)
        else:
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

    # Also do an animation
    plotter = Plotter(coords, energies, "ΔE / au", interval=250)
    plotter.animate()
    plt.show()


def load_results(keys):
    if isinstance(keys, str):
        keys = (keys, )
    results_fn = "results.yaml"
    with open("results.yaml") as handle:
        all_results = yaml.load(handle.read())
    num_cycles = len(all_results)

    results_list = list()
    for key in keys:
        tmp_list = list()
        for res_per_cycle in all_results:
            try:
                tmp_list.append([res[key] for res in res_per_cycle])
            except KeyError:
                print(f"Key '{key}' not present in {results_fn}. Exiting.")
                sys.exit()
        results_list.append(np.array(tmp_list))
    # The length of the second axis correpsonds to the number of images
    num_images = results_list[0].shape[1]
    # Flatten the first axis when we got only a single key
    if len(results_list) == 1:
        results_list = results_list[0]
    print(f"Found path with {num_images} images.")
    print(f"Loaded {num_cycles} cycle(s).")
    return results_list, num_cycles, num_images


def plot_cosgrad():
    keys = ("energy", "forces", "coords")
    (energies, forces, coords), num_cycles, num_images = load_results(keys)
    dummy_atoms = list()

    all_nebs = list()
    all_perp_forces = list()
    for i, per_cycle in enumerate(zip(energies, forces, coords), 1):
        ens, frcs, crds = per_cycle
        images = [Geometry(dummy_atoms, per_image) for per_image in crds]
        neb = NEB(images)
        neb._forces = frcs
        neb.forces = frcs
        neb._energy = ens
        neb.energy = ens
        all_nebs.append(neb)
        pf = neb.perpendicular_forces.reshape(num_images, -1)
        all_perp_forces.append(pf)
    """
    last_neb = all_nebs[-1]
    for img in last_neb.images:
        #print(last_forces.shape)
        max_force = img.forces.max()
        rms_force = np.sqrt(np.mean(np.square(img.forces)))
        print(f"rms(force)={rms_force:.04f}, max(force)={max_force:.04}")
    last_pf = all_perp_forces[-1]
    for per_img in last_pf:
        max_force = per_img.max()
        rms_force = np.sqrt(np.mean(np.square(per_img.flatten())))
        print(f"max(force)={max_force:.06}, rms(force)={rms_force:.06f}")
    max_force = last_pf[1:].max()
    rms_force = np.sqrt(np.mean(np.square(last_pf[1:].flatten())))
    """
    all_max_forces = list()
    all_rms_forces = list()
    rms = lambda arr: np.sqrt(np.mean(np.square(arr)))
    for pf in all_perp_forces:
        max_forces = pf.max(axis=1)
        all_max_forces.append(max_forces)
        rms_forces = np.apply_along_axis(rms, 1, pf)
        all_rms_forces.append(rms_forces)
    all_max_forces = np.array(all_max_forces)
    all_rms_forces = np.array(all_rms_forces)

    cmap = plt.get_cmap("Greys")
    fig, (ax1, ax2) = plt.subplots(2)
    # Using dataframes seems to be the easiest way to include
    # the colormap... Axis.plot() cant be used with cmap.
    max_df = pd.DataFrame(all_max_forces.T)
    rms_df = pd.DataFrame(all_rms_forces.T)
    kwargs = {
        "colormap": cmap,
        "logy": True,
        "legend": False,
        #"ylim": (5e-4, 1e-1),
        "marker": "o",
        "xticks": range(num_images),
        "xlim": (0, num_images-1),
    }
    ax1 = max_df.plot(
            ax=ax1,
            title="max(perpendicular grad.)",
            **kwargs,
    )
    ax1.set_xlabel("Image")
    ax2 = rms_df.plot(
            ax=ax2,
            title="rms(perpendicular grad.)",
            **kwargs,
    )
    ax2.set_xlabel("Image")
    plt.tight_layout()
    plt.show()


def plot_multistate_pes(keys):
    (pes_ens, coords), num_cycles, num_images = load_results(keys)
    pes_ens -= pes_ens.min(axis=(2, 1), keepdims=True)

    plotter = Plotter(coords, pes_ens, "ΔE / au")
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
    coords, num_cycles, num_images = load_results(key)
    ac = list()
    for i, per_cycle in enumerate(coords):
        pc = list()
        for j, per_image in enumerate(per_cycle):
            coords_slice = per_image.reshape(-1, 3)[inds]
            param = func(coords_slice)
            pc.append(param)
        ac.append(pc)
    ac_arr = np.array(ac)

    plotter = Plotter(coords, ac_arr, ylabel)
    plotter.animate()

    df = pd.DataFrame(ac_arr.T)
    cmap = plt.get_cmap("Greys")
    ax = df.plot(
            title=f"Params {inds}",
            colormap=cmap,
            legend=False,
            marker="o",
            xticks=range(num_images),
            xlim=(0, num_images-1),
    )
    ax.set_xlabel("Image")
    ax.set_ylabel(ylabel)
    plt.tight_layout()
    plt.show()


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--until", type=int,
                        help="Only show until cycle [until].")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--saras", action="store_true",
                       help="Plot OpenMolcas state average potential energy "
                            "surfaces over the course of the NEB.")
    group.add_argument("--tddft", action="store_true",
                       help="Plot ORCA TDDFT potential energy surfaces "
                            "over the course of the NEB.")
    group.add_argument("--params", type=int, nargs="+")
    group.add_argument("--cosgrad", action="store_true",
                        help="Plot image gradients along the path.")
    group.add_argument("--energies", action="store_true",
                        help="Plot energies.")

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    if args.energies:
        plot_energies()
    elif args.saras:
        keys = ("sa_energies", "coords")
        plot_multistate_pes(keys)
    elif args.tddft:
        keys = ("tddft_energies", "coords")
        plot_multistate_pes(keys)
    elif args.params:
        plot_params(args.params)
    elif args.cosgrad:
        plot_cosgrad()


if __name__ == "__main__":
    run()
