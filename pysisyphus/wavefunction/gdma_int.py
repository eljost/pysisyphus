from typing import Tuple

import numpy as np

from pysisyphus.numint import MolGrid
from pysisyphus.wavefunction import density_numba, Shells, Wavefunction


def get_prim_data(shells: Shells) -> Tuple[np.ndarray, np.ndarray]:
    nprims = shells.nprims
    # Two integers; angular momentum and cartesian starting index
    Ls_inds = np.empty((nprims, 2), dtype=int)
    # Five doubles; contraction coefficient, orbital exponent and three center coordinates
    prim_data = np.empty((nprims, 5), dtype=float)

    i = 0
    for shell in shells:
        L = shell.L
        cart_index = shell.index
        center = shell.center
        for coeff, exp_ in zip(shell.coeffs, shell.exps):
            prim_data[i] = coeff, exp_, *center
            Ls_inds[i] = L, cart_index
            i += 1
    return Ls_inds, prim_data


def get_diffuse_density_numba(
    wf: Wavefunction, mol_grid: MolGrid, switch: float
) -> np.ndarray:
    # Convert pysisyphus shells to simple arrays that can be passed to numba
    Ls_inds, prim_data = get_prim_data(wf.shells)

    rho_pseudo = np.empty_like(mol_grid.weights)

    # Convert (spherical) density matrix to Cartesian matrix
    P_tot = wf.P_tot
    reorder_c2s = wf.shells.reorder_c2s_coeffs
    P_tot_cart = reorder_c2s.T @ P_tot @ reorder_c2s

    density_numba.eval_prim_density(
        Ls_inds,
        prim_data,
        mol_grid.xyz,
        P_tot_cart,
        switch,
        rho_pseudo,
    )
    return rho_pseudo
