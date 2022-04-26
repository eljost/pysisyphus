import argparse
from pathlib import Path
import sys
import textwrap

import h5py
import matplotlib
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import splrep, splev

from pysisyphus.constants import AU2KJPERMOL, AU2KCALPERMOL, AU2EV
from pysisyphus.config import OUT_DIR_DEFAULT
from pysisyphus.dynamics import Gaussian
from pysisyphus.io import parse_xyz
from pysisyphus.peakdetect import peakdetect
from pysisyphus.wrapper.jmol import render_cdd_cube


CDD_PNG_FNS = "cdd_png_fns"
"""
Default to kJ mol⁻¹. DE_LABEL and get_en_conv will be overwritten when
run() is executed. Both are still defined here, to make the functions
from plot.py importable.
"""
DE_LABEL = r"$\Delta$E / kJ mol⁻¹"


def get_en_conv():
    return AU2KJPERMOL, "kJ mol⁻¹"


def get_force_unit(coord_type):
    force_unit = "$E_h$ Bohr⁻¹"
    if coord_type != "cart":
        force_unit += " (rad)⁻¹"
    return force_unit


def spline_plot_cycles(cart_coords, energies):
    num_cycles = energies.shape[1]

    fig, ax = plt.subplots()
    colors = matplotlib.cm.Greys(np.linspace(0.2, 1, num=num_cycles))
    # for cycle, color in zip(energies, colors):
    for i, (cycle, color) in enumerate(zip(energies, colors)):
        ax.plot(cycle, "o-", color=color)
    ax.set_title("COS image energies")

    kwargs = {
        "ls": ":",
        "color": "darkgrey",
    }
    # Try to spline the last cycle to get an estimate for the spliend HEI
    try:
        last_cycle = energies[-1]
        num_images = last_cycle.size
        spl = splrep(np.arange(num_images), last_cycle)
        # Calculate interpolated values
        x2 = np.linspace(0, num_images, 100)
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
    _, en_unit = get_en_conv()
    ax.set_ylabel(DE_LABEL)

    return fig, ax


def plot_cycle(cart_coords, energies):
    # Plot last_cycle
    fig, ax = plt.subplots()
    last_energies = energies[-1].copy()
    xs = np.arange(len(last_energies))
    ax.plot(xs, last_energies, "o-")
    ax.set_xlabel("Image")
    ax.set_ylabel(DE_LABEL)
    ax.set_title(f"COS image energies, (last) cycle {len(energies)-1}")

    first_image_en = last_energies[0]
    last_image_en = last_energies[-1]
    max_en_ind = np.nanargmax(last_energies)
    max_en = last_energies[max_en_ind]
    nans = np.isnan(last_energies).sum()
    if nans:
        print(
            f"The COS seems to be not fully grown (yet), {nans} energies are missing."
        )
    print(
        "Barrier heights using actual energies (not splined) from "
        f"cycle {energies.shape[0]-1}."
    )
    print(f"\tHighest energy image (HEI) at index {max_en_ind} (0-based)")

    first_barr = max_en - first_image_en
    _, en_unit = get_en_conv()
    print(f"\tBarrier between first image and HEI: {first_barr:.1f} {en_unit}")
    last_barr = max_en - last_image_en
    print(f"\tBarrier between last image and HEI: {last_barr:.1f} {en_unit}")

    return fig, ax


def anim_cos(cart_coords, energies):
    num_cycles = cart_coords.shape[0]

    # Also do an animation
    min_ = np.nanmin(energies)
    max_ = np.nanmax(energies)

    coord_diffs = np.linalg.norm(cart_coords - cart_coords[0][0], axis=2)
    fig, ax = plt.subplots()

    # Initial energies
    lines = ax.plot(coord_diffs[0], energies[0], "o-")
    y_max = max_ - min_
    ax.set_ylim(0, y_max)
    ax.set_xlabel("Coordinate differences / Bohr")
    ax.set_ylabel(DE_LABEL)

    def update_func(i):
        fig.suptitle("Cycle {}".format(i))
        lines[0].set_xdata(coord_diffs[i])
        lines[0].set_ydata(energies[i])

    def animate():
        animation = FuncAnimation(
            fig,
            update_func,
            frames=num_cycles,
            interval=250,
        )
        return animation

    anim = animate()
    return anim, fig, ax


def load_h5(h5_fn, h5_group, datasets=None, attrs=None):
    if datasets is None:
        datasets = list()

    if attrs is None:
        attrs = list()

    with h5py.File(h5_fn, "r") as handle:
        group = handle[h5_group]

        atoms = group.attrs["atoms"]
        cur_cycle = group.attrs["cur_cycle"]
        coord_size = group.attrs["coord_size"]
        num_cycles = cur_cycle + 1

        image_nums = group["image_nums"][:num_cycles].astype(int)
        image_inds = group["image_inds"][:num_cycles].astype(int)

        _datasets = dict()
        for ds in datasets:
            try:
                _datasets[ds] = group[ds][:num_cycles]
            except KeyError:
                print(f"Could not load dataset '{ds}' from HDF5 file.")

        _attrs = dict()
        for a in attrs:
            try:
                _attrs[a] = group.attrs[a]
            except KeyError:
                print(f"Could not load attribute '{a}' from HDF5 file.")

    en_conv, _ = get_en_conv()
    if "energies" in _datasets:
        ens = _datasets["energies"]
        ens -= ens.min()
        ens *= en_conv

    try:
        # We can't use coord_size because coord_type may be != cart and then
        # coord_size gives the number of internals.
        cart_shape = (num_cycles, -1, 3 * len(atoms))
        _datasets["cart_coords"] = _datasets["cart_coords"].reshape(cart_shape)
    except KeyError:
        pass

    try:
        # Here we can use coord_size because forces will always be in the same
        # coordinate system as the actual coordinates.
        _datasets["forces"] = _datasets["forces"].reshape((num_cycles, -1, coord_size))
    except KeyError:
        pass

    def sort_by_image(arr):
        by_image = np.full_like(arr, np.nan)
        for cyc, (img_ind, img_num) in enumerate(zip(image_inds, image_nums)):
            img_ind = img_ind[:img_num]
            by_image[cyc, img_ind] = arr[cyc, :img_num]
        return by_image

    for k, v in _datasets.items():
        _datasets[k] = sort_by_image(v)

    # Also copy requested attributes into dictionary
    _datasets.update(_attrs)

    return _datasets


def plot_cos_energies(h5_fn="optimization.h5", h5_group="opt"):
    results = load_h5(
        h5_fn, h5_group, datasets=("cart_coords", "energies"), attrs=("is_cos",)
    )
    cart_coords = results["cart_coords"]
    energies = results["energies"]

    assert results["is_cos"]

    # Splined last cycle and plot of all cycles
    fig_, ax_ = spline_plot_cycles(
        cart_coords, energies
    )  # lgtm [py/unused-local-variable]
    # Plot last cycle
    fig_last, ax_last = plot_cycle(
        cart_coords, energies
    )  # lgtm [py/unused-local-variable]
    # Plot animation
    anim, fig_anim, ax_anim = anim_cos(
        cart_coords, energies
    )  # lgtm [py/unused-local-variable]

    plt.show()


def plot_cos_forces(h5_fn="optimization.h5", h5_group="opt", last=15):
    results = load_h5(
        h5_fn,
        h5_group,
        datasets=("energies", "forces",),
        attrs=("is_cos", "coord_type", "max_force_thresh", "rms_force_thresh"),
    )
    cycles = len(results["energies"])
    last_cycles = np.arange(cycles)[-last:]
    energies = results["energies"][-last:]
    forces = results["forces"][-last:]
    coord_type = results["coord_type"]

    assert results["is_cos"]

    last_axis = forces.ndim - 1
    max_ = np.nanmax(np.abs(forces), axis=last_axis)
    rms = np.sqrt(np.mean(forces ** 2, axis=last_axis))
    hei_indices = energies.argmax(axis=1)
    force_unit = get_force_unit(coord_type)

    fmt = ".6f"
    print("HEI forces in E_h / a0")
    for i, hei_index in enumerate(hei_indices):
        cycle = last_cycles[i]
        hei_max = max_[i, hei_index]
        hei_rms = rms[i, hei_index]
        print(f"\tCycle {cycle:03d}: max(forces)={hei_max:{fmt}}, rms(forces)={hei_rms:{fmt}}")

    fig, (ax0, ax1) = plt.subplots(sharex=True, nrows=2)

    def plot(ax, data, title):
        num = data.shape[0]
        alphas = np.linspace(0.125, 1, num=num)
        colors = matplotlib.cm.Greys(np.linspace(0, 1, num=num))
        colors[-1] = (1., 0., 0., 1.)  # use red for latest cycle
        for row, color, alpha in zip(data, colors, alphas):
            ax.plot(row, "o-", color=color, alpha=alpha)
            ax.set_ylabel(force_unit)
        ax.set_yscale("log")
        if title:
            ax.set_title(title)

    plot(ax0, max_, "max(perp. forces)")
    plot(ax1, rms, "rms(perp. forces)")
    try:
        ax0.axhline(results["max_force_thresh"], ls="--", c="k", label="Conv. thresh.")
        ax0.legend()
        ax1.axhline(results["rms_force_thresh"], ls="--", c="k", label="Conv. thresh.")
        ax1.legend()
    except KeyError:
        print("Could not find max/rms entries for force threshold on HDF5 file!")
    ax1.set_xlabel("Image")

    plt.tight_layout()
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
        steps.append(i + 1)
    else:
        print("Root flip occured in the last step. Not showing the last step.")

    energies = np.array(energies_)
    roots = np.array(roots_)

    fig, ax = plt.subplots()
    for i, state in enumerate(energies.T):
        ax.plot(steps, state, "o-", label=f"State {i:03d}")
    ax.legend(loc="lower center", ncol=3)
    ax.set_xlabel("Cycle")
    ax.set_ylabel(r"$\Delta$E / eV")
    root_ens = [s[r] for s, r in zip(energies, roots)]
    ax.plot(steps, root_ens, "--k")
    plt.show()


def plot_md(h5_group="run"):
    with h5py.File("md.h5", "r") as handle:
        group = handle[h5_group]

        steps = group["step"][:]
        ens = group["energy_tot"][:]
        ens_conserved = group["energy_conserved"][:]
        Ts = group["T"][:]
        T_avgs = group["T_avg"][:]
        # coords = group["cart_coords"][:]
        # velocities = group["velocity"][:]

        dt = group.attrs["dt"]
        T_target = group.attrs["T_target"]

    _, (ax0, ax1, ax2) = plt.subplots(nrows=3, sharex=True)
    dts = steps * dt
    en_conv, en_unit = get_en_conv()
    ens *= en_conv
    mean = ens.mean()
    ens -= mean
    ax0.plot(dts, ens)
    ax0.axhline(0, ls="--", c="k")
    ax0.set_ylabel(r"$E - \overline{E}$ / " + en_unit)
    ax0.set_title("Energy")

    ens_conserved *= en_conv
    mean_conserved = ens_conserved.mean()
    ens_conserved -= mean_conserved
    ax1.plot(dts, ens_conserved)
    ax1.axhline(0, ls="--", c="k")
    ax1.set_ylabel(r"$E_\\mathrm{cons.} - \overline{E}_\\mathrm{cons.}$ / {en_unit}")
    ax1.set_title("Conserved quantity")

    ax2.plot(dts, Ts, label="Current")
    ax2.plot(dts, T_avgs, ls="--", c="orange", label="Average")
    ax2.axhline(T_target, ls="--", c="k", label="Target")
    ax2.legend()
    ax2.set_title(f"mean(T) = {Ts.mean():.2f} K")
    ax2.set_xlabel(r"$\Delta t$ / fs")
    ax2.set_ylabel("T / K")

    plt.tight_layout()
    plt.show()


def plot_gau(gau_fns, num=50):
    print("Assuming constant Gaussian s & w!")

    assert (
        0 < len(gau_fns) < 3
    ), "Currently, only plotting of 1 or 2 collective variables is possible!"

    gaussians = list()
    centers = list()
    grids = list()
    for gau_fn in gau_fns:
        gau_data = np.loadtxt(gau_fn)
        gau_centers = gau_data[:, 3]
        _, s, w, _ = gau_data[0]
        gau = Gaussian(w=w, s=s)
        print(f"Successfully loaded '{gau_fn}' with w={w:.6f}, s={s:.6f}")
        gaussians.append(gau)
        min_ = gau_centers.min()
        max_ = gau_centers.max()
        diff = abs(max_ - min_)
        centers.append(gau_centers)
        grid = np.linspace(min_, max_, num=num)
        grids.append(grid)
        print(f"\tmin={min_:.4f}, max={max_:.4f}, Δ={diff:.4f}")

    def eval_gaussians(coords):
        value = 0.0
        for i, gau in enumerate(gaussians):
            value += gau.value(coords=coords[i], x0=centers[i])
        return value

    fig, ax = plt.subplots()

    en_conv, en_unit = get_en_conv()
    if len(gau_fns) == 1:
        grid = grids[0]
        ens = -np.array([eval_gaussians(x) for x in grid[:, None]]) * en_conv
        ens -= ens.min()
        ax.plot(grid, ens)
        ax.set_xlabel(f"CV0, {gau_fns[0]}")
        ax.set_ylabel(rf"$\Delta F$ / {en_unit}")
    elif len(gau_fns) == 2:
        grid0, grid1 = grids
        X, Y = np.meshgrid(grid0, grid1)
        xy_flat = np.stack((X.flatten(), Y.flatten()), axis=1)
        ens = (
            -np.array([eval_gaussians(xy) for xy in xy_flat]).reshape(num, num)
            * en_conv
        )
        ens -= ens.min()
        levels = np.linspace(ens.min(), 0.75 * ens.max(), num=15)
        # contour = ax.contour(X, Y, ens, levels=levels)
        _ = ax.contourf(X, Y, ens, levels=levels)
        # plt.clabel(contour, inline=True, fmt="%1.1f", fontsize=10, colors="white", levels=levels)
        ax.set_xlabel(f"CV0, {gau_fns[0]}")
        ax.set_ylabel(f"CV1, {gau_fns[1]}")
    plt.show()


def plot_overlaps(h5, thresh=0.1):
    with h5py.File(h5, "r") as handle:
        overlaps = handle["overlap_matrices"][:]
        roots = handle["roots"][:]
        calculated_roots = handle["calculated_roots"][:]
        ref_cycles = handle["ref_cycles"][:]
        ref_roots = handle["ref_roots"][:]
        try:
            ovlp_type = handle.attrs["ovlp_type"]
            ovlp_with = handle.attrs["ovlp_with"]
        # The old way is handled below. Newer pysis versions store ovlp_type/ovlp_with
        # in attrs.
        except KeyError:
            ovlp_type = handle["ovlp_type"][()].decode()
            ovlp_with = handle["ovlp_with"][()].decode()
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
        print(f"Found rendered {len(cdd_imgs)} CDD images.")

    overlaps[np.abs(overlaps) < thresh] = np.nan
    print(f"Overlap type: {ovlp_type}")
    print(f"Overlap with: {ovlp_with}")
    print(f"Found {len(overlaps)} overlap matrices.")
    print(f"Roots: {roots}")
    print(f"Reference cycles: {ref_cycles}")
    print(f"Reference roots: {ref_roots}")
    print()
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
        ax.imshow(o, vmin=0, vmax=1)
        ax.grid(color="#CCCCCC", linestyle="--", linewidth=1)
        ax.set_xticks(np.arange(n_states, dtype=np.int))
        ax.set_yticks(np.arange(n_states, dtype=np.int))
        # set_ylim is needed, otherwise set_yticks drastically shrinks the plot
        ax.set_ylim(n_states - 0.5, -0.5)
        ax.set_xlabel("new roots")
        ax.set_ylabel("reference roots")
        for (l, k), value in np.ndenumerate(o):
            if np.isnan(value):
                continue
            value_str = f"{abs(value):.2f}"
            ax.text(k, l, value_str, ha="center", va="center")
        j, k = ref_cycles[i], i + 1
        ref_root = ref_roots[i]
        ref_ind = ref_root - 1
        if ovlp_type == "wf":
            ref_ind += 1
        old_root = calculated_roots[i + 1]
        new_root = roots[i + 1]
        ref_overlaps = o[ref_ind]
        argmax = np.nanargmax(ref_overlaps)
        xy = (argmax - 0.5, ref_ind - 0.5)
        highlight = Rectangle(xy, 1, 1, fill=False, color="red", lw="4")
        ax.add_artist(highlight)
        if ax1:
            ax1.imshow(cdd_imgs[i])
        fig.suptitle(
            f"overlap {i:03d}\n"
            f"{ovlp_type} overlap between {j:03d} and {k:03d}\n"
            f"old root: {old_root}, new root: {new_root}"
        )
        fig.canvas.draw()

    draw(0)

    i = 0
    i_backup = i
    i_last = len(overlaps) - 1

    def press(event):
        nonlocal i
        nonlocal i_backup
        if event.key == "left":
            i = max(0, i - 1)
        elif event.key == "right":
            i = min(i_last, i + 1)
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

    plt.tight_layout()
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
    png_stems = [png.stem for png in png_fns if png.exists()]
    print(f"{len(png_stems)} cubes seem already rendered.")

    # Only render cubes that are not yet rendered
    cdd_cubes = [cube for cube in cdd_cubes if Path(cube).stem not in png_stems]
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


def plot_afir(h5_fn="afir.h5", h5_group="afir"):

    h5_fns = (h5_fn, Path(OUT_DIR_DEFAULT) / h5_fn)
    for h5_fn in h5_fns:
        print(f"Trying to open '{h5_fn}' ... ", end="")
        try:
            with h5py.File(h5_fn, "r") as handle:
                group = handle[h5_group]
                cycles = group.attrs["cur_cycle"] + 1
                afir_ens = group["energy"][:cycles]
                true_ens = group["true_energy"][:cycles]
                afir_forces = group["forces"][:cycles]
                true_forces = group["true_forces"][:cycles]
            print("done.")
            break
        except FileNotFoundError:
            print("file not found.")
            continue


    en_conv, en_unit = get_en_conv()
    afir_ens *= en_conv
    afir_ens -= afir_ens.min()
    true_ens *= en_conv
    true_ens -= true_ens.min()
    afir_forces = np.linalg.norm(afir_forces, axis=1)
    true_forces = np.linalg.norm(true_forces, axis=1)

    fig, (en_ax, forces_ax) = plt.subplots(nrows=2, sharex=True)

    style1 = "r--"
    style2 = "g--"
    style3 = "bo-"

    l1 = en_ax.plot(afir_ens, style1, label="AFIR")
    l2 = en_ax.plot(true_ens, style2, label="True")
    en_ax2 = en_ax.twinx()
    l3 = en_ax2.plot(true_ens + afir_ens, style3, label="Sum")
    en_ax2.tick_params(axis="y", labelcolor="blue")

    lines = l1 + l2 + l3
    labels = [l.get_label() for l in lines]
    en_ax.legend(lines, labels, loc=0)

    en_ax.set_title("Energies")
    en_ax.set_ylabel(DE_LABEL)

    forces_ax.set_title("||Forces||")
    l1 = forces_ax.plot(afir_forces, style1, label="AFIR")
    l2 = forces_ax.plot(true_forces, style2, label="True")

    forces_ax2 = forces_ax.twinx()
    l3 = forces_ax2.plot(true_forces + afir_forces, style3, label="Sum")
    forces_ax2.tick_params(axis="y", labelcolor="blue")

    lines = l1 + l2 + l3
    labels = [l.get_label() for l in lines]
    forces_ax.legend(lines, labels, loc=0)
    forces_ax.set_xlabel("Cycle")
    forces_ax.set_ylabel("$E_h$ Bohr⁻¹")

    peak_inds, _ = peakdetect(true_ens, lookahead=2)
    if peak_inds:
        print("Peaks:")
        print(f"\tCycle: Energy / {en_unit}")
        print()
        for at, energy in peak_inds:
            print(f"\t{at}: {energy:.2f}")

    try:
        peak_xs, peak_ys = zip(*peak_inds)
        highest = np.argmax(peak_ys)

        en_ax.scatter(peak_xs, peak_ys, s=100, marker="X", c="k", zorder=10)
        en_ax.scatter(
            peak_xs[highest], peak_ys[highest], s=150, marker="X", c="k", zorder=10
        )
        en_ax.axvline(peak_xs[highest], c="k", ls="--")
        forces_ax.axvline(peak_xs[highest], c="k", ls="--")
    except ValueError as err:
        print("Peak-detection failed!")

    # fig.legend(loc="upper right")
    plt.tight_layout()
    plt.show()


def parse_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5_fn", default="overlap_data.h5")
    parser.add_argument("--h5_group", default="opt", help="HDF5 group to plot.")
    parser.add_argument("--orient", default="")
    parser.add_argument(
        "--kcal", action="store_true", help="Use kcal mol⁻¹ instead of kJ mol⁻¹."
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--cosforces",
        "--cf",
        action="store_true",
        help="Plot image forces along a COS.",
    )
    group.add_argument(
        "--cosens", "--ce", action="store_true", help="Plot COS energies."
    )
    group.add_argument(
        "--all_energies",
        "-a",
        action="store_true",
        help="Plot ground and excited state energies from 'overlap_data.h5'.",
    )
    group.add_argument(
        "--afir",
        action="store_true",
        help="Plot AFIR and true -energies and -forces from an AFIR calculation.",
    )
    group.add_argument("--opt", action="store_true", help="Plot optimization progress.")
    group.add_argument("--irc", action="store_true", help="Plot IRC progress.")
    group.add_argument("--overlaps", "-o", action="store_true")
    group.add_argument("--render_cdds", action="store_true")
    group.add_argument("--h5_list", default=None, help="List groups in HDF5 file.")
    group.add_argument("--md", action="store_true", help="Plot MD.")
    group.add_argument("--gau", nargs="*")
    group.add_argument("--scan", action="store_true")
    group.add_argument(
        "--trj", help="Plot energy values from the comments of a " ".trj file."
    )

    return parser.parse_args(args)


def plot_opt(h5_fn="optimization.h5", h5_group="opt"):
    with h5py.File(h5_fn, "r") as handle:
        try:
            group = handle[h5_group]
        except KeyError:
            groups = list(handle.keys())
            groups_str = "\t" + "\n\t".join(groups)
            print(
                f"Could not find group '{h5_group}'!\nAvailable groups are:\n{groups_str}\n"
                f"Use '--h5_group [group]' to plot a different group."
            )
            if groups:
                group = handle[groups[0]]
                print(f"Using first group '{group}'.")
            else:
                return

        cur_cycle = group.attrs["cur_cycle"]
        is_cos = group.attrs["is_cos"]
        is_converged = group.attrs["is_converged"]
        coord_type = group.attrs["coord_type"]

        ens = group["energies"][:cur_cycle]
        max_forces = group["max_forces"][:cur_cycle]
        rms_forces = group["rms_forces"][:cur_cycle]
        max_force_thresh = group.attrs["max_force_thresh"]
        rms_force_thresh = group.attrs["rms_force_thresh"]

    en_conv, en_unit = get_en_conv()
    ens -= ens.min()
    ens *= en_conv
    if is_cos:
        text = textwrap.wrap(
            "COS optimization detected. Plotting total energy of all images "
            "in every cycle. Results from optimizing growing COS methods can "
            "be plotted but the plots are not really useful as the varying "
            "number of images is not considered.",
            width=80,
        )
        print("\n".join(text))
        ens = ens.sum(axis=1)
    force_unit = get_force_unit(coord_type)

    ax_kwargs = {
        "marker": "o",
    }

    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, sharex=True)

    ax0.plot(ens, **ax_kwargs)
    ax0.set_ylabel(DE_LABEL)

    ax1.plot(max_forces, **ax_kwargs)
    ax1.axhline(max_force_thresh, c="k", ls="--", label="Threshold")
    ax1.set_yscale("log")
    ax1.set_title("max(forces)")
    ax1.set_ylabel(force_unit)
    ax1.legend()

    ax2.plot(rms_forces, **ax_kwargs)
    ax2.axhline(rms_force_thresh, c="k", ls="--", label="Threshold")
    ax2.set_yscale("log")
    ax2.set_title("rms(forces)")
    ax2.set_xlabel("Cycle")
    ax2.set_ylabel(force_unit)
    ax2.legend()

    title = f"{h5_fn}/{h5_group}, converged={is_converged}"
    fig.suptitle(title, y=0.999)

    plt.tight_layout()
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
    with h5py.File(h5, "r") as handle:
        mw_coords = handle["mw_coords"][:]
        energies = handle["energies"][:]
        gradients = handle["gradients"][:]
        rms_grad_thresh = handle["rms_grad_thresh"][()]

        try:
            ts_index = handle["ts_index"][()]
        except KeyError:
            ts_index = None

    sizes = [dataset.shape[0] for dataset in (mw_coords, energies, gradients)]
    size0 = sizes[0]
    assert all([size == size0 for size in sizes])
    print(f"\tFound {size0} IRC points")

    en_conv, en_unit = get_en_conv()
    energies -= energies[0]
    energies *= en_conv

    cds = np.linalg.norm(mw_coords - mw_coords[0], axis=1)
    rms_grads = np.sqrt(np.mean(gradients ** 2, axis=1))
    max_grads = np.abs(gradients).max(axis=1)

    fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, sharex=True)

    plt_kwargs = {
        "linestyle": "-",
        "marker": "o",
    }

    ax0.plot(cds, energies, **plt_kwargs)
    ax0.set_title("energy change")
    ax0.set_ylabel(DE_LABEL)

    ax1.plot(cds, rms_grads, **plt_kwargs)
    ax1.axhline(rms_grad_thresh, linestyle="--", color="k")
    ax1.set_title("rms(gradient)")
    ax1.set_ylabel("$E_h$ Bohr⁻¹")

    ax2.plot(cds, max_grads, **plt_kwargs)
    ax2.set_title("max(gradient)")
    ax2.set_xlabel("IRC / amu$^{\\frac{1}{2}}$ Bohr")
    ax2.set_ylabel("$E_h$ Bohr⁻¹")

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


def plot_scan(h5_fn="scan.h5"):
    with h5py.File(h5_fn, "r") as handle:
        groups = list()
        energies = list()
        for k in handle.keys():
            group = handle[k]
            groups.append(k)
            energies.append(group["energies"][:])
    print(f"Found {len(groups)} groups in '{h5_fn}'.")

    en_conv, _ = get_en_conv
    for group, ens in zip(groups, energies):
        ens -= ens.min()
        ens *= en_conv
        fig, ax = plt.subplots()
        ax.plot(ens, "o-")
        ax.set_xlabel("Scan point")
        ax.set_ylabel(DE_LABEL)
        ax.set_title(group)
    plt.show()


def plot_trj_energies(trj):
    """Parse comments of .xyz/.trj as energies and plot."""
    atoms_coords, comments = parse_xyz(trj, with_comment=True)
    try:
        energies = np.array(comments, dtype=float)
    except ValueError as err:
        print("Could not convert comments to energies!\n")
        raise err
    en_conv, en_unit = get_en_conv()
    energies -= energies.min()
    energies *= en_conv
    fig, ax = plt.subplots()
    ax.plot(energies)
    highlights = [0, energies.argmax(), energies.size - 1]
    highlight_ens = energies[highlights]
    ax.scatter(highlights, highlight_ens, s=50)
    ax.set_xlabel("Step")
    ax.set_ylabel(DE_LABEL)
    ax.set_title(trj)
    plt.show()


def list_h5_groups(h5_fn):
    with h5py.File(h5_fn, "r") as handle:
        groups = list(handle.keys())

    print(f"Found {len(groups)} groups in '{h5_fn}'\n")
    for i, grp in enumerate(groups):
        print(f"\t{i:02d}: {grp}")

    if groups:
        print("\nAvailable groups can be selected by '--h5_group [name]'.")


def run():
    args = parse_args(sys.argv[1:])

    h5_fn = Path(args.h5_fn)

    global get_en_conv

    def get_en_conv():
        if args.kcal:
            conv, label = AU2KCALPERMOL, "kcal mol⁻¹"
        else:
            conv, label = AU2KJPERMOL, "kJ mol⁻¹"
        return conv, label

    global DE_LABEL
    DE_LABEL = rf"$\Delta$E / {get_en_conv()[1]}"

    if args.all_energies or args.overlaps:
        if not h5_fn.exists():
            if (tmp_fn := Path(OUT_DIR_DEFAULT) / h5_fn).exists():
                h5_fn = tmp_fn
                print(f"Loading overlap data from '{h5_fn}'.")

    # Optimization
    if args.h5_list:
        list_h5_groups(args.h5_list)
    if args.opt:
        plot_opt(h5_group=args.h5_group)
    # COS specific
    elif args.cosens:
        plot_cos_energies(h5_group=args.h5_group)
    elif args.cosforces:
        plot_cos_forces(h5_group=args.h5_group)
    # AFIR
    elif args.afir:
        plot_afir()
    # IRC related
    elif args.irc:
        plot_irc()
    # Overlap calculator related
    elif args.all_energies:
        plot_all_energies(h5=h5_fn)
    elif args.overlaps:
        plot_overlaps(h5=h5_fn)
    # MD related
    elif args.md:
        plot_md()
    elif args.gau:
        plot_gau(args.gau)
    elif args.scan:
        plot_scan()
    elif args.render_cdds:
        render_cdds(h5=h5_fn)
    elif args.trj:
        plot_trj_energies(args.trj)


if __name__ == "__main__":
    run()
