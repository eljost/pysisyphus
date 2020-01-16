#!/usr/bin/env python3

from pyscf.tools import cubegen
from pyscf.lib.chkfile import load_mol


from pysisyphus.helpers import geom_from_library
from pysisyphus.calculators.PySCF import PySCF


def get_density_cube(mol, mf, cube_fn="density.cub"):
    cubegen.density(mol, cube_fn, mf.make_rdm1())


def get_esp_cube(mol, mf, cube_fn="esp.cub"):
    cubegen.mep(mol, cube_fn, mf.make_rdm1())


def gen_cubes():
    geom = geom_from_library("benzene.xyz")
    calc = PySCF(basis="sto-3g")
    geom.set_calculator(calc)

    geom.energy

    # mol = load_mol("calculator_000.000.mol.chk")
    # # Density
    # get_density_cube(mol, calc.mf)
    # # ESP, extremly slow, lol
    # # https://github.com/pyscf/pyscf-doc/blob/master/examples/1-advanced/031-MEP.py
    # get_esp_cube(mol, calc.mf)


if __name__ == "__main__":
    gen_cubes()
