from pathlib import Path

import numpy as np

from pysisyphus.config import CONFIG_DIR


NUMINT_DIR = Path(CONFIG_DIR / "numint")
_LEBEDEV_GRIDS = dict()


def get_lebedev_grid(n):
    """Get Cartesian Lebedev-grid with selected number of points.

    Parameters
    ----------
    n
        Number of quadrature points.

    Returns
    -------
    2d array of shape (n, 4) containing the angular Lebedev grid. The first
    three columns contain the Cartesian coordinates of the grid points on the
    unit sphere. The last column contains the weights. They are already multiplied
    by 4π!
    """
    # Try to lookup grid in the dict
    if not n in _LEBEDEV_GRIDS:
        # Weights already include factor of 4π
        grids = np.load(NUMINT_DIR / "lebedev_grids_4pi.npz")
        # Grids are stored with string-keys in the file
        grid = grids[str(n)]
        _LEBEDEV_GRIDS[n] = grid
    # Always return a grid copy
    return _LEBEDEV_GRIDS[n].copy()
