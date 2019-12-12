#!/usr/bin/env python3

import os
from pathlib import Path

import numpy as np

from pysisyphus.helpers import geom_from_xyz_file
from pysisyphus.stocastic.align import matched_rmsd

THIS_DIR = Path(os.path.dirname(os.path.realpath(__file__)))


def test_matched_rmsd():
    geom1 = geom_from_xyz_file(THIS_DIR / "eins.xyz")
    # Calling with the identical geometries should return RMSD of 0.
    min_rmsd, (geom1_matched, geom2_matched) = matched_rmsd(geom1, geom1)
    np.testing.assert_allclose(min_rmsd, 0.0, atol=1e-10)
    np.testing.assert_allclose(geom1_matched.coords, geom2_matched.coords)

    geom2 = geom_from_xyz_file(THIS_DIR / "zwei.xyz")
    min_rmsd, _ = matched_rmsd(geom1, geom2)
    np.testing.assert_allclose(min_rmsd, 0.057049, atol=1e-5)


if __name__ == "__main__":
    test_matched_rmsd()
