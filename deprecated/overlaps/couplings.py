#!/usr/bin/env python3

import itertools as it

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d


from pysisyphus.overlaps.sorter import sort_cut

np.set_printoptions(suppress=True, precision=2)

def interpolate(x, ys, x_fine, **kwargs):
    if x is None:
        x = np.arange(ys.size)
    f = interp1d(x, ys, **kwargs)
    y_fine = f(x_fine)
    return y_fine


def get_state_couplings(ad_ens, dia_ens, ind_a, ind_b, x=None):
    assert ind_a != ind_b
    
    if x is None:
        x = np.arange(ad_ens.shape[0])
    x_fine = np.linspace(x.min(), x.max(), 512)
    interp = lambda ys: interpolate(x=x, ys=ys, x_fine=x_fine)
    a_ad = ad_ens[:,ind_a]
    b_ad = ad_ens[:,ind_b]
    a_ad_fine = interp(a_ad)
    b_ad_fine = interp(b_ad)
    ad_fine = interpolate(x, dia_ens, x_fine, axis=0)

    a_dia = dia_ens[:,ind_a]
    b_dia = dia_ens[:,ind_b]
    a_dia_fine = interp(a_dia)
    b_dia_fine = interp(b_dia)
    dia_fine = interpolate(x, dia_ens, x_fine, axis=0)

    frac = np.abs(
            (a_dia_fine - b_ad_fine) / (a_ad_fine - b_ad_fine)
    )
    # frac_root = np.full_like(frac)
    frac_root = np.sqrt(frac)
    # Cap at 1
    frac_root[frac_root > 1] = 1
    alpha = np.arccos(frac_root)
    # print(frac_root)
    # print(frac_root.max())
    Vd = np.abs(
        (b_ad_fine - a_ad_fine) * np.cos(alpha) * np.sin(alpha)
    )

    grey = "#aaaaaa"

    fig, (ax_ad, ax_dia, ax_Vd) = plt.subplots(nrows=3)
    title = f"Couplings between State {ind_a+1} and {ind_b+1}"
    ax_ad.plot(x_fine, ad_fine, c=grey)
    ax_ad.plot(x_fine, a_ad_fine)
    ax_ad.plot(x_fine, b_ad_fine)
    ax_ad.set_title("adiabatic")

    ax_dia.plot(x_fine, dia_fine, c=grey)
    ax_dia.plot(x_fine, a_dia_fine)
    ax_dia.plot(x_fine, b_dia_fine)
    ax_dia.set_title("diabatic")

    ax_Vd.plot(x_fine, Vd)
    ax_Vd.set_title("couplings")

    fig.suptitle(title)

    # plt.tight_layout()
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])

    plt.show()


def run():
    # ens from parse_tddft.py
    ens = np.loadtxt("all_energies.dat")
    # Drop GS
    ens = ens[:,1:]

    # fig, ax = plt.subplots()
    # ax.plot(ens, "o-")
    # plt.show()

    max_ovlp_inds = np.loadtxt("tden_max_ovlp_inds", dtype=int)
    # print("max_overlap_inds")
    # print(max_ovlp_inds)

    consider_first = 3
    ens_sorted, inds_sorted = sort_cut(ens, max_ovlp_inds,
                                       consider_first=consider_first)

    # fig, ax = plt.subplots()
    # ax.plot(ens_sorted, "o-")
    # plt.show()

    get_state_couplings(ens, ens_sorted, 2, 1)
    #get_state_couplings(ens, ens_sorted, 1, 2)
    get_state_couplings(ens, ens_sorted, 2, 0)
    # get_state_couplings(ens, ens_sorted, 0, 2)


def d2d3():
    adia = np.loadtxt("d2d3_adia.csv")
    dia = np.loadtxt("d2d3_dia.csv")
    x = adia[:,0]
    adia = adia[:,1:]
    dia = dia[:,1:]
    get_state_couplings(adia, dia, 0, 1, x=x)
    # get_state_couplings(adia, dia, 1, 0, x=x)
    get_state_couplings(adia, dia, 1, 2, x=x)
    # get_state_couplings(adia, dia, 2, 1, x=x)


def couplings(states):
    energies = np.loadtxt("all_energies.dat")
    # Drop GS
    energies = energies[:,1:]
    max_ovlp_inds = np.loadtxt("tden_max_ovlp_inds", dtype=int)
    combs = it.combinations(states, 2)
    # coupling_mat = [get_state_couplings
    consider_first = max(states)
    energies_sorted, inds_sorted = sort_cut(energies, max_ovlp_inds,
                                            consider_first=consider_first)
    # fig, ax = plt.subplots()
    # ax.plot(ens_sorted, "o-")
    # plt.show()
    state_couplings = {
        (ind_a, ind_b): get_state_couplings(
                            energies, energies_sorted, ind_a, ind_b
                        ) for ind_a, ind_b in combs
    }

if __name__ == "__main__":
    run()
    # d2d3()
