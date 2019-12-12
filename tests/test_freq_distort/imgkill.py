#!/usr/bin/env python3

from collections import namedtuple
import sys

import matplotlib.pyplot as plt
import numpy as np

# from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.calculators.XTB import XTB
from pysisyphus.helpers import geom_from_xyz_file, eigval_to_wavenumber
from pysisyphus.TablePrinter import TablePrinter
from pysisyphus.xyzloader import coords_to_trj


np.set_printoptions(suppress=True, precision=4)

NMResult = namedtuple("NMResult", "energy gradient hessian w v nus neg_inds neg_nus")


def do_analysis(geom, nu_thresh=-5.):
    # Do calculations
    hess = geom.hessian
    gradient = geom.gradient
    energy = geom.energy

    # Calculate hessian
    hessian = geom.hessian
    mw_hess = geom.mw_hessian

    proj_hess = geom.eckart_projection(mw_hess)
    w, v = np.linalg.eigh(proj_hess)
    nus = eigval_to_wavenumber(w)
    neg_inds = nus <= nu_thresh

    grad_norm = np.linalg.norm(gradient)
    print(f"@norm(grad)={grad_norm:.6f}")
    # Overlaps between gradient and imaginary modes. We use the absolute
    # values, as the sign of the eigenvectors is ambiguous.
    overlaps = np.einsum("ij,i->j", v[:,neg_inds], gradient)

    print(f"Found {neg_inds.sum()} imaginary frequencies < {nu_thresh} cm⁻¹.")
    tp = TablePrinter("# nu/cm⁻¹ <grad|qi>".split(), "int float float".split())
    tp.print_header()
    for i, (nu, ovlp) in enumerate(zip(nus[neg_inds], overlaps)):
        tp.print_row((i, nu, ovlp))
    res = NMResult(
            energy=energy,
            gradient=gradient,
            hessian=hess,
            w=w,
            v=v,
            nus=nus,
            neg_inds=neg_inds,
            neg_nus=nus[neg_inds],
    )
    return res


def relax_along_mode(geom, mode, step_size=.1, max_steps=25, nu_thresh=-5.):
    org_coords = geom.coords.copy()
    hessian = geom.hessian
    mw_hess = geom.mw_hessian


    mw_grad = geom.mw_gradient
    overlap = mode @ mw_grad
    if overlap > 0.:
        mode *= -1
    lengths = np.arange(max_steps) * step_size
    coords_list = list()
    M = np.sqrt(geom.masses_rep)
    results = list()
    mw_coords = geom.mw_coords.copy()
    results = list()
    used_lengths = list()
    energies = list()
    converged = False
    coords = list()
    print("Running energy calculations ", end="")
    for l in lengths:
        print(".", end="")
        new_coords = (mw_coords + l*mode) / M
        coords.append(new_coords.copy())
        geom.coords = new_coords
        energy = geom.energy
        sys.stdout.flush()
        energies.append(energy)
        used_lengths.append(l)
        coords_list.append(new_coords)
        try:
            if energy > energies[-2]:
                print(" energy increased!")
                break
        except IndexError:
            pass
    else:
        print(" max cycles reached!")

    energies = np.array(energies)
    coords = np.array(coords)
    used_lengths = np.array(used_lengths)
    return energies, coords, used_lengths, converged


def run_relaxations(geom, mode_ind, nu_thresh=-5.):
    res = do_analysis(geom)
    print(f"Distorting along mode #{mode_ind}")
    org_mode = res.v[:,mode_ind]
    mode = org_mode.copy()
    for i in range(3):
        mode = res.v[:,mode_ind]
        _ = relax_along_mode(geom, mode, nu_thresh=nu_thresh)
        energies, coords, used_lengths, converged = _

        # geom.coords = coords[-2]

        energies -= energies.min()
        energies *= 2625.50

        fig, ax = plt.subplots()
        ax.plot(used_lengths, energies, "o--")
        ax.set_ylabel("$\Delta$E / kJ mol⁻¹")
        fig.suptitle(f"Cycle {i:02d}")

        plt.tight_layout()
        plt.show()
        fig.savefig(f"imgkill_{i:02d}.pdf")
        res = do_analysis(geom)
        overlaps = abs(np.einsum("ij,i->j", res.v, org_mode))
        max_ovlp_ind = overlaps.argmax()
        print( "Original imargy mode has highest overlap with mode "
              f"{max_ovlp_ind} ({overlaps[max_ovlp_ind]:.2%})"
        )
        eigenvalue = res.w[max_ovlp_ind]
        nu = eigval_to_wavenumber(eigenvalue)
        print(f"@Wavenumber of mode {mode_ind:02d}: {nu:.1f} cm⁻¹")
        if nu > nu_thresh:
            break
        print()


def run():
    geom = geom_from_xyz_file("03_00_water_addition_ts.xyz")
    geom.set_calculator(XTB(pal=4))
    run_relaxations(geom, 1)


if __name__ == "__main__":
    run()
