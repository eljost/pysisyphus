import ctypes as ct
from pathlib import Path
from typing import Callable, Tuple

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


def get_fortran_eval_prim_density_func() -> Callable:
    this_dir = Path(__file__)
    flib = ct.CDLL(this_dir.with_name("./gdma_int_fortran.so"))
    eval_prim_density = flib.eval_prim_density
    eval_prim_density.argtypes = [
        ct.c_int32,  # nprims
        ct.c_int32,  # npoints
        ct.c_int32,  # nbfs
        ct.POINTER(ct.c_int32),  # Ls_inds
        ct.POINTER(ct.c_double),  # primdata
        ct.POINTER(ct.c_double),  # coords3d
        ct.POINTER(ct.c_double),  # P
        ct.c_double,  # switch
        ct.POINTER(ct.c_double),  # rho
    ]
    return eval_prim_density


def i32ptr(arr: np.ndarray):
    return arr.ctypes.data_as(ct.POINTER(ct.c_int32))


def dblptr(arr: np.ndarray):
    return arr.ctypes.data_as(ct.POINTER(ct.c_double))


def get_diffuse_density_fortran(
    wf: Wavefunction, mol_grid: MolGrid, switch: float
) -> np.ndarray:
    eval_prim_density = get_fortran_eval_prim_density_func()

    # Convert pysisyphus shells to simple arrays that can be passed to Fortran.
    Ls_inds, primdata = get_prim_data(wf.shells)

    rho_pseudo = np.empty_like(mol_grid.weights)
    coords3d = mol_grid.xyz
    P_tot = wf.P_tot
    nprims = len(primdata)
    npoints = len(mol_grid.xyz)

    # Convert (spherical) density matrix to Cartesian matrix
    reorder_c2s = wf.shells.reorder_c2s_coeffs
    P_tot_cart = reorder_c2s.T @ P_tot @ reorder_c2s
    nbfs = P_tot_cart.shape[0]

    # Make 2d arrays F-contiguous
    Ls_inds_f = np.asfortranarray(Ls_inds.astype(np.int32).T)
    primdata_f = np.asfortranarray(primdata.T)
    coords3d_f = np.asfortranarray(coords3d.T)
    P_tot_cart_f = np.asfortranarray(P_tot_cart)

    eval_prim_density(
        ct.c_int32(nprims),
        ct.c_int32(npoints),
        ct.c_int32(nbfs),
        #
        i32ptr(Ls_inds_f),
        dblptr(primdata_f),
        dblptr(coords3d_f),
        dblptr(P_tot_cart_f),
        ct.c_double(switch),
        dblptr(rho_pseudo),
    )
    return rho_pseudo
