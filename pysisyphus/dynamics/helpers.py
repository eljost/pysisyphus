#!/usr/bin/env python3

import numpy as np

from pysisyphus.xyzloader import make_trj_str
from pysisyphus.constants import BOHR2ANG


def dump_coords(atoms, coords, trj_fn):
    coords = np.array(coords)
    coords = coords.reshape(-1, len(atoms), 3) * BOHR2ANG
    trj_str = make_trj_str(atoms, coords)

    with open(trj_fn, "w") as handle:
        handle.write(trj_str)
