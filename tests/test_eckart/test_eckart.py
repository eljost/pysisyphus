#!/usr/bin/env python3

import numpy as np

from pysisyphus.helpers import geom_from_library
from pysisyphus.calculators.Gaussian16 import Gaussian16


def eigvals2nus(vals):
    """https://chemistry.stackexchange.com/questions/74639"""
    C = 137  # Speed of light in atomic units
    TO_NU = 1 / (2*np.pi*C*5.29177e-9)
    el_mass = 5.4858e-4 # electron mass in atomic units
    nus = np.sign(vals) * np.sqrt(np.abs(vals)*el_mass) * TO_NU
    for i, n in enumerate(nus):
        print(f"{i:02d}: {n:.4f}")
    return nus


def test_eckart():
    geom = geom_from_library("benzene_pm6_opt.xyz")
    calc = Gaussian16(route="PM6", pal=4)
    geom.set_calculator(calc)
    H = geom.hessian
    np.savetxt("bz_hess", H)

    # H = np.loadtxt("bz_hess")
    # geom.hessian = H
    # vals, vecsT = np.linalg.eigh(geom.mw_hessian)
    # eigvals2nus(vals)
    # print()

    H_mw_proj = geom.eckart_projection(geom.mw_hessian)
    valsp, vecsTp = np.linalg.eigh(H_mw_proj)
    eigvals2nus(valsp)
    assert (np.abs(valsp[:6]) < 1e-16).all()


if __name__ == "__main__":
    test_eckart()
