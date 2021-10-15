#!/usr/bin/env python3

import argparse
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d


from pysisyphus.constants import AU2EV
from pysisyphus.helpers import index_array_from_overlaps

np.set_printoptions(suppress=True, precision=4)


def sort_cut(cut, max_overlap_inds, consider_first=None):
    # For N-1 overlaps we have N points
    assert max_overlap_inds.shape[0] == cut.shape[0]-1
    if consider_first is None:
        consider_first = cut.shape[1]
    cut_sorted = cut.copy()
    cur_inds = np.arange(cut.shape[1])
    # Start from index 1, as we keep the first row
    all_inds = [cur_inds, ]
    for i, inds in enumerate(max_overlap_inds, 1):
        jumps = [(j, s) for j, s in enumerate(inds) if s != j]

        unique_inds = np.unique(inds[:consider_first])
        # if jumps:
            # print(f"from {i-1:02d} to {i:02d}", jumps)

        new_inds = cur_inds.copy()
        if unique_inds.size != consider_first:
            print(f"Step {i:02d}, indices are non-unique! Not skipping!")
            print(f"\t ({i-1} -> {i}): {inds}")
        else:
            for j, s in jumps:
                new_inds[j] = cur_inds[s]
        
        cur_inds = new_inds
        cut_sorted[i:,cur_inds] = cut[i:]
        all_inds.append(cur_inds)
    all_inds = np.array(all_inds)
    return cut_sorted, all_inds


def interpolate(x, ys, x_fine, **kwargs):
    if x is None:
        x = np.arange(ys.size)
    f = interp1d(x, ys, **kwargs)
    y_fine = f(x_fine)
    return y_fine


def get_state_couplings(ad_ens, dia_ens, ind_a, ind_b, x=None):
    assert ind_a != ind_b

    if x is None:
        x = np.arange(ad_ens.shape[0])
    x_fine = np.linspace(x.min(), x.max(), 512)
    interp = lambda ys: interpolate(x=x, ys=ys, x_fine=x_fine)
    a_ad = ad_ens[:, ind_a]
    b_ad = ad_ens[:, ind_b]
    a_ad_fine = interp(a_ad)
    b_ad_fine = interp(b_ad)
    ad_fine = interpolate(x, dia_ens, x_fine, axis=0)

    a_dia = dia_ens[:, ind_a]
    b_dia = dia_ens[:, ind_b]
    a_dia_fine = interp(a_dia)
    b_dia_fine = interp(b_dia)
    dia_fine = interpolate(x, dia_ens, x_fine, axis=0)

    frac = np.abs((a_dia_fine - b_ad_fine) / (a_ad_fine - b_ad_fine))
    # frac_root = np.full_like(frac)
    frac_root = np.sqrt(frac)
    # Cap at 1
    frac_root[frac_root > 1] = 1.0
    alpha = np.arccos(frac_root)
    Vd = np.abs((b_ad_fine - a_ad_fine) * np.cos(alpha) * np.sin(alpha))

    grey = "#aaaaaa"

    fig, (ax_ad, ax_dia, ax_Vd) = plt.subplots(nrows=3, sharex=True)
    title = f"Couplings between State {ind_a+1} and {ind_b+1}"
    ax_ad.plot(x_fine, ad_fine, c=grey)
    ax_ad.plot(x_fine, a_ad_fine)
    ax_ad.plot(x_fine, b_ad_fine)
    ax_ad.set_title("adiabatic")

    ax_dia.plot(x_fine, dia_fine, c=grey)
    ax_dia.plot(x_fine, a_dia_fine)
    ax_dia.plot(x_fine, b_dia_fine)
    ax_dia.set_title("diabatic")

    ax_Vd.plot(x_fine, Vd)
    ax_Vd.set_title("couplings")

    fig.suptitle(title)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

    return x_fine, Vd


def from_overlap_data(h5_fn, ind_a, ind_b):
    with h5py.File(h5_fn, "r") as handle:
        all_energies = handle["all_energies"][:]
        root_flips = handle["root_flips"][:]
        overlap_mats = handle["overlap_matrices"][:]
        ovlp_type = handle.attrs["ovlp_type"]

    print("Root flips:")
    for i, rf in enumerate(root_flips):
        print("Step", i, rf)

    all_energies -= all_energies.min()
    all_energies *= AU2EV
    ad_ens = all_energies

    index_arr = [index_array_from_overlaps(ovlp_mat) for ovlp_mat in overlap_mats]
    index_arr = np.array(index_arr, dtype=int)
    if ovlp_type != "wf":
        steps, states = index_arr.shape
        ext_shape = (steps, states + 1)
        _index_arr = np.zeros(ext_shape, dtype=int)
        _index_arr[:, 1:] = index_arr + 1
        index_arr = _index_arr

    dia_ens, _ = sort_cut(all_energies, index_arr)
    ad_ens = ad_ens[~root_flips]
    dia_ens = dia_ens[~root_flips]

    return get_state_couplings(ad_ens, dia_ens, ind_a=ind_a, ind_b=ind_b)


def parse_args(args):
    parser = argparse.ArgumentParser()

    parser.add_argument("states", type=int, nargs=2)
    parser.add_argument("--h5_fn", default="overlap_data.h5")

    return parser.parse_args(args)


def run():
    args = parse_args(sys.argv[1:])

    ind_a, ind_b = args.states
    h5_fn = args.h5_fn
    from_overlap_data(h5_fn, ind_a, ind_b)


if __name__ == "__main__":
    run()
