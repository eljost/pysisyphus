#!/usr/bin/env python3

import re

import numpy as np

from pysisyphus.helpers import geom_from_library
from pysisyphus.constants import BOHR2ANG, ANG2BOHR


def print_gaussian_ints(geom):
    def get_line(base, index, definition, value):
        def_str = "(" + ",".join([str(d+1) for d in definition]) + ")"
        return f" ! {base}{index:<5d} {base}{def_str:<21} {value:.4f}"
    coord_dict = dict()
    bonds = 1
    angles = 1
    dihedrals = 1
    string_list = []
    for pc in geom.internal._prim_internals:
        if len(pc.inds) == 2:
            base = "R"
            val = pc.val*BOHR2ANG
            index = bonds
            bonds += 1
        elif len(pc.inds) == 3:
            base = "A"
            val = np.rad2deg(pc.val)
            index = angles
            angles += 1
        elif len(pc.inds) == 4:
            base = "D"
            val = np.rad2deg(pc.val)
            index = dihedrals
            dihedrals += 1
        coord_dict[tuple(sorted(pc.inds+1))] = val
        string_list.append(get_line(base, index, pc.inds, val))
    return string_list, coord_dict


def read_gaussian_ints(text):
    coord_dict = dict()
    for line in text.strip().split("\n"):
        _, name, definition, value, *rest = line.split()
        value = float(value)
        inds = re.match("[RAD]\(([\d,]+)\)", definition)[1].split(",")
        #inds_set = frozenset(tuple([int(i) for i in inds]))
        #print(inds_set)
        inds_tpl = tuple(sorted([int(i) for i in inds]))
        #print(definition, value)
        #if len(inds) == 2:
        #    value *= ANG2BOHR
        #else:
        #    value = np.deg2rad(value)
        #coord_dict[inds_set] = value
        if inds_tpl in coord_dict:
            print(f"its already there! {inds_tpl}")
        coord_dict[inds_tpl] = value
    return coord_dict


def compare_to_gaussian(ints1, ints2):
    for inds1, val1 in ints1.items():
        try:
            val2 = ints2[inds1]
            diff = val1-val2
            print(f"{val1:>12.4f}{val2:>12.4f}{diff:>12.4f}")
        except KeyError:
            print(f"{inds1} not in ints2!")
    keyset1 = set(ints1.keys())
    keyset2 = set(ints2.keys())
    print("keys that are only in ints1", keyset1-keyset2)
    print("keys that are only in ints2", keyset2-keyset1)


if __name__ == "__main__":
    fn = "azetidine_gaussian_internals"
    with open(fn) as handle:
        text = handle.read()
    cd = read_gaussian_ints(text)
