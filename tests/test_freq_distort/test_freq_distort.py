#!/usr/bin/env python3

from collections import namedtuple
from pathlib import Path
import sys

import cloudpickle
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from pysisyphus.calculators.ORCA import ORCA
from pysisyphus.calculators.XTB import XTB
from pysisyphus.helpers import geom_from_xyz_file, eigval_to_wavenumber, shake_coords
from pysisyphus.irc.ModeKill import ModeKill
from pysisyphus.TablePrinter import TablePrinter
from pysisyphus.xyzloader import coords_to_trj


np.set_printoptions(suppress=True, precision=4)


def from_orca(path="."):
    calc = ORCA("")
    path = Path(path)
    results = calc.parse_hessian(path)
    hess = results["hessian"]
    np.save("hess", hess)
    # hess = np.load("hess.npy")
    return hess


def run():
    geom = geom_from_xyz_file("shaked.xyz")
    calc = XTB(pal=4)
    geom.set_calculator(calc)
    kill_modes(geom)
    # tmp(geom, 1)


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
    print(f"norm(grad)={grad_norm:.6f}")
    print()
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

def tmp(geom, mode_ind, nu_thresh=-5.):
    org_coords = geom.coords.copy()
    hessian = geom.hessian
    mw_hess = geom.mw_hessian
    start_result = do_analysis(geom, nu_thresh)

    print(f"Distorting along mode #{mode_ind}")
    mode = start_result.v[:,mode_ind]

    # lengths = np.linspace(-2, 2, 41)
    lengths = np.linspace(-1.5, 1.5, 31)
    print(lengths)
    coords_list = list()
    M = np.sqrt(geom.masses_rep)
    results = list()
    mw_coords = geom.mw_coords.copy()
    for l in lengths:
        print(f"l={l:02f}")
        new_coords = (mw_coords + l*mode) / M
        geom.coords = new_coords
        # res = do_analysis(geom)
        # sys.stdout.flush()
        # results.append(res)
        
        coords_list.append(new_coords)
    coords_to_trj("displaced.trj", geom.atoms, coords_list)
    # with open("res", "wb") as handle: cloudpickle.dump(results, handle)
    with open("res", "rb") as handle: results=cloudpickle.load(handle)

    energies = np.array([r.energy for r in results])
    energies -= energies.min()
    energies *= 2625.50
    grad_norms = [np.linalg.norm(r.gradient) for r in results]
    neg_nus = np.array([r.nus[:2] for r in results])

    fig, (ax, ax2, ax3) = plt.subplots(nrows=3, sharex=True)
    ax.plot(lengths, energies, "o--")
    ax.set_title("Energies")
    ax.set_ylabel("$\Delta$E / kJ mol⁻¹")
    ax2.plot(lengths, grad_norms, "o--")
    ax2.set_ylabel("Hartree Bohr⁻¹")
    ax2.set_title("norm(gradient)")

    width = .025
    for i, modes in enumerate(neg_nus.T):
        ax3.bar(lengths+i*width, modes, width=width)
    ax3.set_xlabel("l")
    ax3.set_ylabel(r"$\nu$ / cm⁻¹")

    plt.tight_layout()
    plt.show()
    fig.savefig("imgkill.pdf")


def kill_modes(geom):

    do_analysis(geom)
    # kill_inds = (0, 1)
    # mk = ModeKill(geom, kill_inds=kill_inds, max_cycles=5)
    kill_inds = (1, )
    mk = ModeKill(geom, kill_inds=kill_inds)
    mk.run()
    do_analysis(geom)


def freq_distort(geom, mode_ind):
    print(geom)
    hess = geom.hessian
    mw_hess = geom.mass_weigh_hessian(hess)
    # w, v = np.linalg.eigh(mw_hess)
    # nus = eigval_to_wavenumber(w)
    # neg_inds = w < 0
    # print(w[neg_inds])
    # print(nus[neg_inds])

    # shaked = geom.copy()
    # shaked.coords = shake_coords(shaked.coords, scale=0.05, seed=25032018)
    # shaked.set_calculator(XTB(pal=4, base_name="shake"))
    # shess = shaked.hessian
    # sw, sv = np.linalg.eigh(shess)
    # with open("shaked.xyz", "w") as handle:
        # handle.write(shaked.as_xyz())

    print("Eckart projected")
    proj_hess = geom.eckart_projection(mw_hess)
    w, v = np.linalg.eigh(proj_hess)
    nus = eigval_to_wavenumber(w)
    neg_inds = w < 0
    print(w[neg_inds])
    print(nus[neg_inds])

    imag_mode = v[:,mode_ind]
    print(imag_mode)
    lengths = np.linspace(-1, 1, 21)
    print(lengths)
    coords_list = list()
    M = np.sqrt(geom.masses_rep)
    for l in lengths:
        new_coords = (geom.mw_coords + l*imag_mode) / M
        coords_list.append(new_coords)
    coords_to_trj("displaced.trj", geom.atoms, coords_list)

    energies = list()
    for cs in coords_list:
        geom.coords = cs
        energy = geom.energy
        energies.append(energy)
        pass
    print(energies)

    energies = np.array(energies)
    np.save("energies", energies)
    # energies -= energies.min()
    # energies *= 2625.25

    def func_harmonic(x, a, b, c):
        return a + b*(x+c)**2

    fromto = 4
    slice_ = slice(fromto, -fromto+1)
    ydata = energies[slice_]
    xdata = np.arange(energies.size)[slice_]
    popt, pcov = curve_fit(func_harmonic, xdata, ydata)
    print("popt")
    print(popt)
    print("pcov")
    print(pcov)

    fig, ax = plt.subplots()
    ax.plot(energies, "o-", label="True")
    ax.plot(xdata, func_harmonic(xdata, *popt), "--", label="Harmonic")
    ax.legend()
    plt.show()

    pass


def test_freq_distort():
    geom = geom_from_xyz_file("final_geometry.xyz")

    print(geom)
    # calc = ORCA("")
    # print(calc)
    # cwd = Path(".")
    # results = calc.parse_hessian(cwd)
    # print(results)
    # hess = results["hessian"]
    # np.save("hess", hess)
    # hess = np.load("hess.npy")

    hess = np.loadtxt("calculated_final_cart_hessian")
    # print("Non-mass-weighted")
    # w, v = np.linalg.eigh(hess)
    # nus = eigval_to_wavenumber(w)
    # neg_inds = w < 0
    # print(w[neg_inds])
    # print(nus[neg_inds])

    # print("Mass-weighted")
    mw_hess = geom.mass_weigh_hessian(hess)
    # w, v = np.linalg.eigh(mw_hess)
    # nus = eigval_to_wavenumber(w)
    # neg_inds = w < 0
    # print(w[neg_inds])
    # print(nus[neg_inds])

    # shaked = geom.copy()
    # shaked.coords = shake_coords(shaked.coords, scale=0.05, seed=25032018)
    # shaked.set_calculator(XTB(pal=4, base_name="shake"))
    # shess = shaked.hessian
    # sw, sv = np.linalg.eigh(shess)
    # with open("shaked.xyz", "w") as handle:
        # handle.write(shaked.as_xyz())

    print("Eckart projected")
    proj_hess = geom.eckart_projection(mw_hess)
    w, v = np.linalg.eigh(proj_hess)
    nus = eigval_to_wavenumber(w)
    neg_inds = w < 0
    print(w[neg_inds])
    print(nus[neg_inds])

    imag_mode = v[:,0]
    print(imag_mode)
    lengths = np.linspace(-1, 1, 21)
    print(lengths)
    coords_list = list()
    M = np.sqrt(geom.masses_rep)
    for l in lengths:
        new_coords = (geom.mw_coords + l*imag_mode) / M
        coords_list.append(new_coords)
    coords_to_trj("displaced.trj", geom.atoms, coords_list)

    # calc = XTB(pal=4)
    # geom.set_calculator(calc)
    # energies = list()
    # for cs in coords_list:
        # geom.coords = cs
        # energy = geom.energy
        # energies.append(energy)
        # pass
    # print(energies)

    # energies = np.array(energies)
    # np.save("energies", energies)
    energies = np.load("energies.npy")
    energies -= energies.min()
    energies *= 2625.25

    def func_harmonic(x, a, b, c):
        return a + b*(x+c)**2

    fromto = 4
    slice_ = slice(fromto, -fromto+1)
    ydata = energies[slice_]
    xdata = np.arange(energies.size)[slice_]
    popt, pcov = curve_fit(func_harmonic, xdata, ydata)
    print("popt")
    print(popt)
    print("pcov")
    print(pcov)

    fig, ax = plt.subplots()
    ax.plot(energies, "o-", label="True")
    ax.plot(xdata, func_harmonic(xdata, *popt), "--", label="Harmonic")
    ax.legend()
    plt.show()

    pass


if __name__ == "__main__":
    # test_freq_distort()
    # run()
    geom = geom_from_xyz_file("03_00_water_addition_ts.xyz")
    geom.set_calculator(XTB(pal=4))
    # freq_distort(geom, 1)
    tmp(geom, 1)
