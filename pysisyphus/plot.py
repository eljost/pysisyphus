#!/usr/bin/env python3

import argparse
import os
import sys

import h5py
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
import pandas as pd
from scipy.interpolate import splrep, splev
import yaml

from pysisyphus.constants import AU2KJPERMOL, BOHR2ANG, AU2EV
from pysisyphus.cos.NEB import NEB
from pysisyphus.Geometry import Geometry
from pysisyphus.peakdetect import peakdetect


class Plotter:
    def __init__(self, coords, data, ylabel, interval=750, save=None,
                 legend=None):
        self.coords = coords
        self.data = data
        self.ylabel = ylabel
        self.interval = interval
        self.save = save
        self.legend = legend

        # First image of the first cycle
        self.anchor = self.coords[0][0]
        self.cycles = len(self.data)
        self.pause = True

        self.fig, self.ax = plt.subplots()
        self.fig.canvas.mpl_connect('key_press_event', self.on_keypress)

        self.ax.set_xlabel("Path length / Bohr")
        y_min = self.data.min()
        y_max = self.data.max()
        self.ax.set_ylim(y_min, y_max)
        self.ax.set_ylabel(self.ylabel)

        self.coord_diffs = self.get_coord_diffs(self.coords)

        if self.data.ndim == 2:
            self.update_func = self.update_plot
        elif self.data.ndim == 3:
            self.update_func = self.update_plot2

    def get_coord_diffs(self, coords, normalize=False):
        coord_diffs = list()
        for per_cycle in coords:
            tmp_list = [0, ]
            for i in range(len(per_cycle)-1):
                diff = np.linalg.norm(per_cycle[i+1]-per_cycle[i])
                tmp_list.append(diff)
            tmp_list = np.cumsum(tmp_list)
            offset = np.linalg.norm(self.anchor-per_cycle[0])
            tmp_list += offset
            if normalize:
                tmp_list /= tmp_list.max()
            coord_diffs.append(tmp_list)
        return np.array(coord_diffs)

    def update_plot(self, i):
        """Use this when only 1 state is present."""
        self.fig.suptitle("Cycle {}".format(i))
        self.lines[0].set_xdata(self.coord_diffs[i])
        self.lines[0].set_ydata(self.data[i])
        if self.save:
            self.save_png(i)

    def update_plot2(self, i):
        """Use this when several states are present."""
        self.fig.suptitle("Cycle {}".format(i))
        for j, line in enumerate(self.lines):
            line.set_ydata(self.data[i][:, j])
        if self.save:
            self.save_png(i)

    def save_png(self, frame):
        frame_fn = f"step{frame}.png"
        if not os.path.exists(frame_fn):
            self.fig.savefig(frame_fn)

    def animate(self):
        self.lines = self.ax.plot(self.coord_diffs[0], self.data[0], "o-")
        if self.legend:
            self.ax.legend(self.lines, self.legend)
        self.animation = matplotlib.animation.FuncAnimation(
                                                self.fig,
                                                self.update_func,
                                                frames=self.cycles,
                                                interval=self.interval)
        if self.save:
            self.animation.save("animation.gif", writer='imagemagick', fps=5)
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

    if isinstance(num_images, list):
        print("Please use --aneb instead of --energies")
        return

    df = pd.DataFrame(energies)
    cmap = plt.get_cmap("Greys")
    df = df.transpose()
    df -= df.values.min()
    df *= AU2KJPERMOL

    # Static plot of path with equally spaced images
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
    plotter = Plotter(coords, energies, "ΔE / au", interval=250, save=False)
    plotter.animate()
    plt.show()


def plot_aneb():
    keys = ("energy", "coords")
    (energies, coords), num_cycles, num_images = load_results(keys)

    # Use coordinates of the first image in the first cycle as
    # anchor for all following cycles.
    first_coords = coords[0][0]

    coord_diffs = list()
    min_ = 0
    for en, c in zip(energies, coords):
        cd = np.linalg.norm(c - first_coords, axis=1)
        min_ = min(0, min(en))
        coord_diffs.append(cd)

    energies_ = list()
    for en in energies:
        en = np.array(en)
        en -= min_
        en *= 2625.499638
        energies_.append(en)

    fig, ax = plt.subplots()
    # Initial energies
    lines = ax.plot(coord_diffs[0], energies_[0], "o-")

    ax.set_xlabel("Coordinate differences / Bohr")
    ax.set_ylabel("$\Delta$J / kJ $\cdot$ mol$^{-1}$")

    def update_func(i):
        fig.suptitle("Cycle {}".format(i))
        lines[0].set_xdata(coord_diffs[i])
        lines[0].set_ydata(energies_[i])

    def animate():
        animation = matplotlib.animation.FuncAnimation(
                                            fig,
                                            update_func,
                                            frames=num_cycles,
                                            interval=250,
        )
        return animation
    anim = animate()
    plt.show()


def load_results(keys):
    if isinstance(keys, str):
        keys = (keys, )
    image_results_fn = "image_results.yaml"
    print(f"Reading {image_results_fn}")
    with open(image_results_fn) as handle:
        all_results = yaml.load(handle.read(), Loader=yaml.Loader)
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
    # Determine the number of images. If we have the same number of images
    # set num_images to this number. Otherwise return a list containing
    # the number of images.
    num_images = np.array([len(cycle) for cycle in results_list[0]])
    if all(num_images[0] == num_images):
        num_images = num_images[0]
        print(f"Found path with {num_images} images.")
    # Flatten the first axis when we got only a single key
    if len(results_list) == 1:
        results_list = results_list[0]
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

    # Calculate norms of true force
    # Shape (cycles, images, coords)
    force_norms = np.linalg.norm(forces, axis=2)
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
    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    # Using dataframes seems to be the easiest way to include
    # the colormap... Axis.plot() cant be used with cmap.
    max_df = pd.DataFrame(all_max_forces.T)
    rms_df = pd.DataFrame(all_rms_forces.T)
    norm_df = pd.DataFrame(force_norms.T)

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

    ax2 = rms_df.plot(
            ax=ax2,
            title="rms(perpendicular grad.)",
            **kwargs,
    )

    ax3 =norm_df.plot(
            ax=ax3,
            title="norm(true grad.)",
            **kwargs,
    )
    ax3.set_xlabel("Image")

    plt.tight_layout()
    plt.show()


def plot_multistate_pes(keys):
    (pes_ens, coords), num_cycles, num_images = load_results(keys)
    pes_ens -= pes_ens.min(axis=(2, 1), keepdims=True)
    pes_ens *= 27.211396

    plotter = Plotter(coords, pes_ens, "ΔE / eV")
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

    def get_dihedral(coords_slice):
        raise Exception("Not implemented yet!")

    type_dict = {
        2: ("bond length / pm", get_bond_length),
        3: ("angle / °", get_angle),
        4: ("dihedral / °", get_dihedral)
    }
    inds_list = [[int(i) for i in i_.split()] for i_ in inds.split(",")]
    ylabels, funcs = zip(*[type_dict[len(inds)] for inds in inds_list])
    assert all([len(inds_list[i]) == len(inds_list[i+1])
                for i in range(len(inds_list)-1)]), "Can only display " \
            "multiple coordinates of the same type (bond, angle or " \
            "dihedral."
    # Just use the first label because they all have to be the same
    ylabel = ylabels[0]

    key = "coords"
    # only allow same type of coordinate if multiple coordinates are given?
    coords, num_cycles, num_images = load_results(key)

    # Coordinates for all images for all cycles
    ac = list()
    for i, per_cycle in enumerate(coords):
        # Coordinates for all images per cycle
        pc = list()
        for j, per_image in enumerate(per_cycle):
            # Coordinates per ind for all images
            pi = list()
            for inds, func in zip(inds_list, funcs):
                coords_slice = per_image.reshape(-1, 3)[inds]
                param = func(coords_slice)
                pi.append(param)
            pc.append(pi)
        ac.append(pc)

    ac_arr = np.array(ac)

    # Construct legend list
    legend = ["-".join([str(i) for i in inds]) for inds in inds_list]
    plotter = Plotter(coords, ac_arr, ylabel, legend=legend)
    plotter.animate()

    #df = pd.DataFrame(ac_arr)
    #cmap = plt.get_cmap("Greys")
    #ax = df.plot(
    #        title=f"Params {inds}",
    #        colormap=cmap,
    #        legend=False,
    #        marker="o",
    #        xticks=range(num_images),
    #        xlim=(0, num_images-1),
    #)
    #ax.set_xlabel("Image")
    #ax.set_ylabel(ylabel)
    #plt.tight_layout()
    plt.show()


def plot_all_energies():
    dump_fn = "overlap_data.h5"
    with h5py.File(dump_fn) as handle:
        energies = handle["all_energies"][:]
        roots = handle["roots"][:]
        flips = handle["root_flips"][:]
    print(f"Found a total of {len(roots)} steps.")
    print(f"{flips} root flips occured.")

    energies -= energies.min()
    energies *= AU2EV

    # Don't plot steps where flips occured
    # energies = np.concatenate((energies[0][None,:], energies[1:,:][~flips]), axis=0)
    energies_ = list()
    roots_ = list()
    steps = list()
    for i, root_flip in enumerate(flips):
        if root_flip:
            print(f"Root flip occured between {i} and {i+1}.")
            continue
        print(f"Using step {i}")
        energies_.append(energies[i])
        roots_.append(roots[i])
        steps.append(i)
    # Don't append last step if a root flip occured there.
    if not flips[-1]:
        energies_.append(energies[-1])
        roots_.append(roots[-1])
        steps.append(i+1)
    else:
        print("Root flip occured in the last step. Not showing the last step.")

    energies = np.array(energies_)
    roots = np.array(roots_)

    fig, ax = plt.subplots()
    for i, state in enumerate(energies.T):
        ax.plot(steps, state, "o-", label=f"State {i:03d}")
    ax.legend()
    ax.set_xlabel("Step")
    ax.set_ylabel("$\Delta E / eV$")
    root_ens = [s[r] for s, r in zip(energies, roots)]
    ax.plot(steps, root_ens, "--k")
    plt.show()


def plot_overlaps(thresh=.1):
    dump_fn = "overlap_data.h5"
    with h5py.File(dump_fn) as handle:
        overlaps = handle["overlap_matrices"][:]
        ovlp_type = handle["ovlp_type"][()].decode()
        ovlp_with = handle["ovlp_with"][()].decode()
        roots = handle["roots"][:]
    overlaps[np.abs(overlaps) < thresh] = np.nan
    print(f"Found {len(overlaps)} overlap matrices.")
    print(f"Roots: {roots}")

    fig, ax = plt.subplots()

    n_states = overlaps[0].shape[0]
    def between(i):
        if ovlp_with == "first":
            return (0, i+1)
        elif ovlp_with == "previous":
            return (i, i+1)
        else:
            raise Exception("I didn't expect that ;)")

    def draw(i):
        fig.clf()
        ax = fig.add_subplot("111")
        o = overlaps[i]
        im = ax.imshow(o, vmin=0, vmax=1)
        # fig.colorbar(im)
        ax.grid(color="#CCCCCC", linestyle='--', linewidth=1)
        ax.set_xticks(np.arange(n_states, dtype=np.int))
        ax.set_yticks(np.arange(n_states, dtype=np.int))
        ax.set_xlabel("new states")
        ax.set_ylabel("old states")
        for (l,k), value in np.ndenumerate(o):
            if np.isnan(value):
                continue
            value_str = f"{abs(value):.2f}"
            ax.text(k, l, value_str, ha='center', va='center')
        i, j = between(i)
        root_i = roots[i]
        root_j = roots[j]
        fig.suptitle(f"{ovlp_type} overlap between {i:03d} and {j:03d}\n"
                     f"old root: {root_i}, new root: {root_j}")
        fig.canvas.draw()
    draw(0)

    i = 0
    def press(event):
        nonlocal i
        if event.key == "left":
            i = max(0, i-1)
        elif event.key == "right":
            i = min(len(overlaps)-1, i+1)
        else:
            return
        draw(i)
    fig.canvas.mpl_connect("key_press_event", press)
    plt.show()


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--first", type=int,
                        help="Only consider the first [first] cycles.")
    parser.add_argument("--last", type=int,
                        help="Only consider the last [last] cycles.")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--saras", action="store_true",
                       help="Plot OpenMolcas state average potential energy "
                            "surfaces over the course of the NEB.")
    group.add_argument("--tddft", action="store_true",
                       help="Plot ORCA TDDFT potential energy surfaces "
                            "over the course of the NEB.")
    group.add_argument("--params",
                       help="Follow internal coordinates over the course of "
                            "the NEB. All atom indices have to be 0-based. "
                            "Use two indices for a bond, three indices for "
                            "an angle and four indices for a dihedral. "
                            "The indices for different coordinates have to "
                            "be separated by ','.")
    group.add_argument("--cosgrad", "--cg", action="store_true",
                        help="Plot image gradients along the path.")
    group.add_argument("--energies", "-e", action="store_true",
                        help="Plot energies.")
    group.add_argument("--aneb", action="store_true",
                        help="Plot Adaptive NEB.")
    group.add_argument("--all_energies", "-a", action="store_true",
        help="Plot ground and excited state energies from 'overlap_data.h5'."
    )
    group.add_argument("--overlaps", "-o", action="store_true")

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
    elif args.aneb:
        plot_aneb()
    elif args.all_energies:
        plot_all_energies()
    elif args.overlaps:
        plot_overlaps()


if __name__ == "__main__":
    run()
