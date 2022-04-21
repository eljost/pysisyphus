import numpy as np
import pytest

from pysisyphus.calculators.AnaPot import AnaPot
from pysisyphus.calculators.PySCF import PySCF
from pysisyphus.linalg import finite_difference_hessian
from pysisyphus.helpers import geom_loader
from pysisyphus.testing import using


def assert_hessians(geom):
    hessian = geom.hessian

    def grad_func(coords):
        return -geom.get_energy_and_forces_at(coords)["forces"]

    fd_hessian = finite_difference_hessian(geom.coords, grad_func, acc=4)
    np.testing.assert_allclose(hessian, fd_hessian, atol=1.6e-4)


def test_fd_hessian():
    geom = AnaPot().get_saddles(i=0)
    assert_hessians(geom)


@using("pyscf")
def test_water_fd_hessian():
    geom = geom_loader("lib:h2o.xyz")
    calc = PySCF(basis="321g", pal=2)
    geom.set_calculator(calc)
    assert_hessians(geom)
