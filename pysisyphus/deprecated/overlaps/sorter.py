#!/usr/bin/env python3

import logging

import numpy as np


def sort_cut(cut, max_overlap_inds, consider_first=None):
    # For N-1 overlaps we have N points
    assert max_overlap_inds.shape[0] == cut.shape[0]-1
    if not consider_first:
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


def load_ascii_or_bin(arr_fn):
    try:
        arr = np.loadtxt(arr_fn)
        return arr
    except UnicodeDecodeError as err:
        msg = f"{arr_fn} seems to be binary. Trying to load it."
        print(msg)
        # logging.exception(f"{arr_fn} seems to be binary. Trying to load it.")

    arr = np.load(arr_fn, )
    return arr


def sort_by_overlaps(energies_fn, max_ovlp_inds_fn, consider_first=None):
    energies = load_ascii_or_bin(energies_fn)
    max_ovlp_inds = load_ascii_or_bin(max_ovlp_inds_fn).astype(int)
    print("Maximum overlap inds")
    print(max_ovlp_inds)
    # Drop GS
    energies = energies[:,1:]
    assert max_ovlp_inds.shape[0] == energies.shape[0]-1
    energies_sorted, inds_sorted = sort_cut(energies, max_ovlp_inds,
                                            consider_first=consider_first)
    ens_sorted_fn = "energies_sorted"
    np.savetxt(f"{ens_sorted_fn}.dat", energies_sorted)
    np.save(ens_sorted_fn, energies_sorted)
    inds_sorted_fn = "all_inds_sorted"
    np.savetxt(f"{inds_sorted_fn}.dat", inds_sorted, fmt="%d")
    np.save(inds_sorted_fn, inds_sorted)
    print()
    print("Indices, sorted")
    for i, row in enumerate(inds_sorted):
        print(f"Step {i:02d}", row)
