#!/usr/bin/env python3

import argparse
import os
from pathlib import Path
import sys

import h5py
import matplotlib
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import splrep, splev
import yaml

from pysisyphus.constants import AU2KJPERMOL, BOHR2ANG, AU2EV
from pysisyphus.cos.NEB import NEB
from pysisyphus.Geometry import Geometry
from pysisyphus.peakdetect import peakdetect
from pysisyphus.wrapper.jmol import render_cdd_cube


CDD_PNG_FNS = "cdd_png_fns"


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
        self.animation = FuncAnimation(
                            self.fig,
                            self.update_func,
                            frames=self.cycles,
                            interval=self.interval
        )
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
    keys = ("energy", "cart_coords")
    (energies, coords), num_cycles, num_images = load_results(keys)

    if isinstance(num_images, list):
        print("Please use --aneb instead of --energies")
        return

    lengths = np.array([len(e) for e in energies])
    equal_lengths = lengths == lengths[-1]
    # Hack to support growing string calculations
    energies = np.array([e for e, l in zip(energies, equal_lengths) if l])
    coords = np.array([c for c, l in zip(coords, equal_lengths) if l])
    num_cycles, num_images = energies.shape

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
    keys = ("energy", "cart_coords")
    (energies, coords), num_cycles, num_images = load_results(keys)

    # Use coordinates of the first image in the first cycle as
    # anchor for all following cycles.
    first_coords = coords[0][0]

    coord_diffs = list()
    min_ = 0
    max_ = max(energies[0])
    for en, c in zip(energies, coords):
        cd = np.linalg.norm(c - first_coords, axis=1)
        min_ = min(min_, min(en))
        max_ = max(max_, max(en))
        coord_diffs.append(cd)

    energies_ = list()
    au2kJmol = 2625.499638
    for en in energies:
        en = np.array(en)
        en -= min_
        en *= au2kJmol
        energies_.append(en)

    fig, ax = plt.subplots()
    # Initial energies
    lines = ax.plot(coord_diffs[0], energies_[0], "o-")
    y_max = (max_ - min_) * au2kJmol
    ax.set_ylim(0, y_max)

    ax.set_xlabel("Coordinate differences / Bohr")
    ax.set_ylabel("$\Delta$J / kJ $\cdot$ mol$^{-1}$")

    def update_func(i):
        fig.suptitle("Cycle {}".format(i))
        lines[0].set_xdata(coord_diffs[i])
        lines[0].set_ydata(energies_[i])

    def animate():
        animation = FuncAnimation(
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
                print(f"Key '{key}' not present in {image_results_fn}. Exiting.")
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
    dummy_atoms = ["H" for _ in coords[0][0]]

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


def plot_all_energies(h5):
    with h5py.File(h5) as handle:
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
    for i, root_flip in enumerate(flips[:-1]):
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
    ax.legend(loc="lower center", ncol=3)
    ax.set_xlabel("Step")
    ax.set_ylabel("$\Delta E / eV$")
    root_ens = [s[r] for s, r in zip(energies, roots)]
    ax.plot(steps, root_ens, "--k")
    plt.show()


def plot_bare_energies(h5):
    with h5py.File(h5) as handle:
        energies = handle["all_energies"][:]
    print(f"Found a total of {len(energies)} steps.")

    energies -= energies.min()
    energies *= AU2EV
    steps = np.arange(len(energies))

    fig, ax = plt.subplots()
    for i, state in enumerate(energies.T):
        ax.plot(steps, state, "o-", label=f"State {i:03d}")
    ax.legend(loc="lower center", ncol=3)
    ax.set_xlabel("Step")
    ax.set_ylabel("$\Delta E / eV$")
    plt.show()


def plot_overlaps(h5, thresh=.1):
    with h5py.File(h5) as handle:
        overlaps = handle["overlap_matrices"][:]
        ovlp_type = handle["ovlp_type"][()].decode()
        ovlp_with = handle["ovlp_with"][()].decode()
        roots = handle["roots"][:]
        calculated_roots = handle["calculated_roots"][:]
        ref_cycles = handle["ref_cycles"][:]
        ref_roots = handle["ref_roots"][:]
        try:
            cdd_img_fns = handle["cdd_imgs"][:]
        except KeyError:
            print(f"Couldn't find image data in '{h5}'.")
            try:
                with open(CDD_PNG_FNS) as handle:
                    cdd_img_fns = handle.read().split()
                print(f"Found image data in '{CDD_PNG_FNS}'")
            except FileNotFoundError:
                cdd_img_fns = None
    cdd_imgs = None
    if cdd_img_fns is not None:
        try:
            cdd_imgs = [mpimg.imread(fn) for fn in cdd_img_fns]
        except FileNotFoundError:
            png_paths = [Path(fn.decode()).name for fn in cdd_img_fns]
            cdd_imgs = [mpimg.imread(fn) for fn in png_paths]

    overlaps[np.abs(overlaps) < thresh] = np.nan
    print(f"Found {len(overlaps)} overlap matrices.")
    print(f"Roots: {roots}")
    print(f"Reference cycles: {ref_cycles}")
    print(f"Reference roots: {ref_roots}")

    print("Key-bindings:")
    print("i: switch between current and first cycle.")
    print("e: switch between current and last cycle.")

    fig, ax = plt.subplots()

    n_states = overlaps[0].shape[0]

    def draw(i):
        fig.clf()
        if cdd_imgs is not None:
            ax = fig.add_subplot(121)
            ax1 = fig.add_subplot(122)
        else:
            ax = fig.add_subplot(111)
            ax1 = None
        o = np.abs(overlaps[i])
        im = ax.imshow(o, vmin=0, vmax=1)
        # fig.colorbar(im)
        ax.grid(color="#CCCCCC", linestyle='--', linewidth=1)
        ax.set_xticks(np.arange(n_states, dtype=np.int))
        ax.set_yticks(np.arange(n_states, dtype=np.int))
        ax.set_xlabel("new states")
        ax.set_ylabel("reference states")
        for (l,k), value in np.ndenumerate(o):
            if np.isnan(value):
                continue
            value_str = f"{abs(value):.2f}"
            ax.text(k, l, value_str, ha='center', va='center')
        j, k = ref_cycles[i], i+1
        ref_root = ref_roots[i]
        ref_ind = ref_root - 1
        old_root = calculated_roots[i+1]
        new_root = roots[i+1]
        ref_overlaps = o[ref_ind]
        if ovlp_type == "wf":
            ref_ind += 1
        argmax = np.nanargmax(ref_overlaps)
        xy = (argmax-0.5, ref_ind-0.5)
        highlight = Rectangle(xy, 1, 1,
                              fill=False, color="red", lw="4")
        ax.add_artist(highlight)
        if ax1:
            ax1.imshow(cdd_imgs[i])
        fig.suptitle(f"overlap {i:03d}\n"
                     f"{ovlp_type} overlap between {j:03d} and {k:03d}\n"
                     f"old root: {old_root}, new root: {new_root}")
        fig.canvas.draw()
    draw(0)

    i = 0
    i_backup = i
    i_last = len(overlaps)-1
    def press(event):
        nonlocal i
        nonlocal i_backup
        if event.key == "left":
            i = max(0, i-1)
        elif event.key == "right":
            i = min(i_last, i+1)
        # Switch between current and first cycle
        elif event.key == "i":
            if i == 0:
                # Restore previous cycle
                i = i_backup
            else:
                # Save current i and jump to the first cycle/image
                i_backup = i
                i = 0
        # Switch between current and last cycle
        elif event.key == "e":
            if i == i_last:
                # Restore previous cycle
                i = i_backup
            else:
                # Save current i and jump to the first cycle/image
                i_backup = i
                i = i_last
        else:
            return
        draw(i)
    fig.canvas.mpl_connect("key_press_event", press)
    plt.show()


def render_cdds(h5):
    with h5py.File(h5) as handle:
        cdd_cubes = handle["cdd_cubes"][:].astype(str)
        orient = handle["orient"][()].decode()
    cdd_cubes = [Path(cub) for cub in cdd_cubes]
    print(f"Found {len(cdd_cubes)} CDD cube filenames in {h5}")
    # Check if cubes exist
    non_existant_cubes = [cub for cub in cdd_cubes if not cub.exists()]
    existing_cubes = [str(cub) for cub in set(cdd_cubes) - set(non_existant_cubes)]
    if any(non_existant_cubes):
        print("Couldn't find cubes:")
        print("\n".join(["\t" + str(cub) for cub in non_existant_cubes]))
        print("Dropping full path and looking only for cube names.")
        cub_names = [cub.name for cub in non_existant_cubes]
        _ = [cub for cub in cub_names if Path(cub).exists()]
        existing_cubes = existing_cubes + _
        cdd_cubes = existing_cubes

    # Create list of all final PNG filenames
    png_fns = [Path(cube).with_suffix(".png") for cube in cdd_cubes]
    # Check which cubes are already rendered
    png_stems = [png.stem for png in png_fns
                 if png.exists()]
    print(f"{len(png_stems)} cubes seem already rendered.")

    # Only render cubes that are not yet rendered
    cdd_cubes = [cube for cube in cdd_cubes
                 if Path(cube).stem not in png_stems]
    print(f"Rendering {len(cdd_cubes)} CDD cubes.")

    for i, cube in enumerate(cdd_cubes):
        print(f"Rendering cube {i+1:03d}/{len(cdd_cubes):03d}")
        _ = render_cdd_cube(cube, orient=orient)
    joined = "\n".join([str(fn) for fn in png_fns])
    with open(CDD_PNG_FNS, "w") as handle:
        handle.write(joined)
    print("Rendered PNGs:")
    print(joined)
    print(f"Wrote list of rendered PNGs to '{CDD_PNG_FNS}'")


def plot_afir():
    with open("image_results.yaml") as handle:
        res = yaml.load(handle.read(), Loader=yaml.loader.Loader)

    afir_ens = [_["energy"] for _ in res]
    true_ens = [_["true_energy"] for _ in res]
    afir_ens = np.array(afir_ens) * AU2KJPERMOL
    afir_ens -= afir_ens.min()
    true_ens = np.array(true_ens) * AU2KJPERMOL
    true_ens -= true_ens.min()

    afir_forces = np.linalg.norm([_["forces"] for _ in res], axis=1)
    true_forces = np.linalg.norm([_["true_forces"] for _ in res], axis=1)
    afir_forces = np.array(afir_forces)
    true_forces = np.array(true_forces)


    fig, (en_ax, forces_ax) = plt.subplots(nrows=2, sharex=True)

    style1 = "r--"
    style2 = "g--"
    style3 = "bo-"

    l1 = en_ax.plot(afir_ens, style1, label="AFIR")
    l2 = en_ax.plot(true_ens, style2, label="True")
    en_ax2 = en_ax.twinx()
    l3 = en_ax2.plot(true_ens+afir_ens, style3, label="Sum")
    en_ax2.tick_params(axis="y", labelcolor="blue")

    lines = l1 + l2 + l3
    labels = [l.get_label() for l in lines]
    en_ax.legend(lines, labels, loc=0)

    en_ax.set_title("Energies")
    en_ax.set_ylabel("$\Delta$E kJ / mol")

    forces_ax.set_title("||Forces||")
    l1 = forces_ax.plot(afir_forces, style1, label="AFIR")
    l2 = forces_ax.plot(true_forces, style2, label="True")

    forces_ax2 = forces_ax.twinx()
    l3 = forces_ax2.plot(true_forces + afir_forces, style3, label="Sum")
    forces_ax2.tick_params(axis="y", labelcolor="blue")

    lines = l1 + l2 + l3
    labels = [l.get_label() for l in lines]
    forces_ax.legend(lines, labels, loc=0)

    peak_inds, _ = peakdetect(true_ens, lookahead=2)
    print(f"Peaks: {peak_inds}")
    try:
        peak_xs, peak_ys = zip(*peak_inds)
        highest = np.argmax(peak_ys)

        en_ax.scatter(peak_xs, peak_ys, s=100, marker="X", c="k", zorder=10)
        en_ax.scatter(peak_xs[highest], peak_ys[highest],
                    s=150, marker="X", c="k", zorder=10)
        en_ax.axvline(peak_xs[highest], c="k", ls="--")
        forces_ax.axvline(peak_xs[highest], c="k", ls="--")
    except ValueError as err:
        print("Peak-detection failed!")


    # fig.legend(loc="upper right")
    plt.tight_layout()
    plt.show()


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--first", type=int,
                        help="Only consider the first [first] cycles.")
    parser.add_argument("--last", type=int,
                        help="Only consider the last [last] cycles.")
    parser.add_argument("--h5", default="overlap_data.h5")
    parser.add_argument("--orient", default="")

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
    group.add_argument("--bare_energies", "-b", action="store_true",
        help="Plot ground and excited state energies from 'overlap_data.h5'."
    )
    group.add_argument("--afir", action="store_true",
        help="Plot AFIR and true -energies and -forces from an AFIR calculation."
    )
    group.add_argument("--opt", action="store_true",
        help="Plot optimization progress."
    )
    group.add_argument("--irc", action="store_true",
        help="Plot IRC progress."
    )
    group.add_argument("--overlaps", "-o", action="store_true")
    group.add_argument("--render_cdds", action="store_true")

    return parser.parse_args(args)


def plot_opt():
    def load(fn):
        print(f"Reading {fn}")
        with open(fn) as handle:
            res = yaml.load(handle.read(), Loader=yaml.Loader)
        energies = np.array(res["energies"])
        energies -= energies.min()
        energies *= AU2KJPERMOL
        return energies
    ens = load("optimizer_results.yaml")

    fig, ax = plt.subplots()

    ax.plot(ens, "o-", label="Cartesian")
    ax.set_xlabel("Step")
    ax.set_ylabel("$\Delta E$ / kJ mol⁻¹")
    ax.legend()
    fig.tight_layout()
    plt.show()


def plot_irc():
    cwd = Path(".")
    h5s = cwd.glob("*irc_data.h5")
    for h5 in h5s:
        type_ = h5.name.split("_")[0]
        title = f"{type_.capitalize()} IRC data"
        _ = plot_irc_h5(h5, title)
    plt.show()


def plot_irc_h5(h5, title=None):
    print(f"Reading IRC data {h5}")
    with h5py.File(h5) as handle:
        mw_coords = handle["mw_coords"][:]
        energies = handle["energies"][:]
        gradients = handle["gradients"][:]
        rms_grad_thresh = handle["rms_grad_thresh"][()]
        try:
            ts_index = handle["ts_index"][()]
        except KeyError:
            ts_index = None

    energies -= energies[0]
    energies *= AU2KJPERMOL

    cds = np.linalg.norm(mw_coords - mw_coords[0], axis=1)
    rms_grads = np.sqrt(np.mean(gradients**2, axis=1))
    max_grads = np.abs(gradients).max(axis=1)

    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, sharex=True)

    plt_kwargs = {
        "linestyle": "-",
        "marker": "o",
    }

    ax0.plot(cds, energies, **plt_kwargs)
    ax0.set_title("energy change")
    ax0.set_ylabel("kJ mol⁻¹")

    ax1.plot(cds, rms_grads, **plt_kwargs)
    ax1.axhline(rms_grad_thresh, linestyle="--", color="k")
    ax1.set_title("rms(gradient)")
    ax1.set_ylabel("Hartree / bohr")

    ax2.plot(cds, max_grads, **plt_kwargs)
    ax2.set_title("max(gradient)")
    ax2.set_xlabel("IRC / amu$^{\\frac{1}{2}}$ bohr")
    ax2.set_ylabel("Hartree / bohr")

    if ts_index:
        x = cds[ts_index]
        for ax, arr in ((ax0, energies), (ax1, rms_grads), (ax2, max_grads)):
            xy = (x, arr[ts_index])
            ax.annotate("TS", xy, fontsize=12, fontweight="bold")

    if title:
        fig.suptitle(title)
    else:
        fig.tight_layout()

    return fig, (ax0, ax1, ax2)


def run():
    args = parse_args(sys.argv[1:])

    h5 = args.h5

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
        plot_all_energies(h5=h5)
    elif args.overlaps:
        plot_overlaps(h5=h5)
    elif args.render_cdds:
        render_cdds(h5=h5)
    elif args.bare_energies:
        plot_bare_energies(h5=h5)
    elif args.afir:
        plot_afir()
    elif args.opt:
        plot_opt()
    elif args.irc:
        plot_irc()


if __name__ == "__main__":
    run()
