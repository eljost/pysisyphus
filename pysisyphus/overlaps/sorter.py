#!/usr/bin/env python3

import numpy as np


def sort_cut(cut, max_overlap_inds, consider_first=None):
    # For N-1 overlaps we have N points
    assert max_overlap_inds.shape[0] == cut.shape[0]-1 
    if not consider_first:
        consider_first = cut.shape[1]
    cut_sorted = cut.copy()
    cur_inds = np.arange(cut.shape[1])
    print("cur_inds")
    # Start from index 1, as we keep the first row
    all_inds = [cur_inds, ]
    for i, inds in enumerate(max_overlap_inds, 1):
        jumps = [(j, s) for j, s in enumerate(inds) if s != j]
        print(f"from {i-1:02d} to {i:02d}", jumps)

        unique_inds = np.unique(inds[:consider_first])

        if i== 9:
            import pdb; pdb.set_trace()
        new_inds = cur_inds.copy()
        if unique_inds.size != consider_first:
            print("inds are non-unique! not skipping!")
            print(f"\t ({i-1} -> {i}): {inds}")
        else:
            for j, s in jumps:
                new_inds[j] = cur_inds[s]
        
        cur_inds = new_inds
        cut_sorted[i:] = cut[i:,cur_inds]
        print(f"{i:02d}: {cur_inds}")
        all_inds.append(cur_inds)
    all_inds = np.array(all_inds)
    return cut_sorted, all_inds
