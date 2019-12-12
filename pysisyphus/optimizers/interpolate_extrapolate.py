#!/usr/bin/env python3

import numpy as np

from pysisyphus.optimizers.gdiis import gdiis, gediis
from pysisyphus.optimizers.line_search2 import poly_line_search


def interpolate_extrapolate(coords, energies, forces, steps,
                            ref_step=None, err_vecs=None, max_vecs=10,
                            gediis_thresh=1e-2, gdiis_thresh=2.5e-3):

    can_gediis = np.sqrt(np.mean(forces[-1]**2)) < gediis_thresh
    can_diis = (ref_step is not None) and (np.sqrt(np.mean(ref_step**2)) < gdiis_thresh)

    # GDIIS check
    if can_diis and (err_vecs is not None) and (ref_step is not None):
        diis_result = gdiis(err_vecs, coords, forces, ref_step, max_vecs)
    # GEDIIS check
    elif can_gediis:
        diis_result = gediis(coords, energies, forces)
    else:
        diis_result = None

    interpol_step = None
    interpol_forces = None
    interpol_energy = None

    if not diis_result or ((diis_result.name == "GDIIS") and (diis_result.N == 2)):
        cur_energy = energies[-1]
        prev_energy = energies[-2]
        cur_grad = -forces[-1]
        prev_grad = -forces[-2]
        prev_step = steps[-1]
        interpol_step, interpol_gradient, interpol_energy = poly_line_search(
                                                                cur_energy, prev_energy,
                                                                cur_grad, prev_grad,
                                                                prev_step, coords,
        )
        if ((interpol_step is not None)
            and (interpol_gradient is not None)
            and (interpol_energy is not None)):
            # prev_coords = coords[-2]
            # new_coords = prev_coords + step
            # geom.coords = new_coords
            interpol_forces = -interpol_gradient
            # interpol_step = step
    elif diis_result:
        interpol_forces = diis_result.forces
        # geom.coords = diis_result.coords
        # Set interpol step
        import pdb; pdb.set_trace()
    return interpol_step, interpol_forces, interpol_energy
